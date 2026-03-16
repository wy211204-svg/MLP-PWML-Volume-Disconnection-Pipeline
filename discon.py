import os
import glob
import shutil
import pandas as pd
import ants
import SimpleITK as sitk
from collections import defaultdict
from utils import run_command_in_shell, get_fsl_volume, PipelineStepEmptyError, logger

class DisconnectionPipeline:
    def __init__(self, config, dirs):
        self.config = config
        self.dirs = dirs

    def run(self):
        logger.info("Disconnection: 开始分析...")
        registered_lesions = self._step1_registration()
        probtrackx_base_dir = self._step2_probtrackx(registered_lesions)
        final_patient_intersections = self._intermediate_processing_chain(probtrackx_base_dir)

        all_patient_ids = [p['id'] for p in self.config['patients_to_process']]
        discon_df = self._step6_calculate_results(final_patient_intersections, all_patient_ids)

        logger.info("--- Disconnection 分析完成 ---")
        if 'Discon. %' in discon_df.columns:
            discon_df['Discon_Score_Numeric'] = pd.to_numeric(
                discon_df['Discon. %'].astype(str).str.replace(r' \(.*\)', '', regex=True),
                errors='coerce'
            ).fillna(0) / 100.0
        return discon_df

    def _step1_registration(self):
        logger.info("Disconnection-Step 1: 运行ANTS配准...")
        registered_lesion_info = []
        controls_fa_t1_dir = self.config['disconnection_inputs']['controls_t1_dir']
        control_files = sorted(
            glob.glob(os.path.join(controls_fa_t1_dir, '*.nii.gz')) + glob.glob(os.path.join(controls_fa_t1_dir, '*.nii')))

        for p_info in self.config['patients_to_process']:
            patient_id, patient_fa_path, lesion_path = p_info['id'], p_info['fa'], p_info['lesion']
            logger.info(f"--- 开始处理患者 {patient_id} 的配准 ---")

            patient_fa_ants = ants.image_read(patient_fa_path)
            patient_lesion_ants = ants.image_read(lesion_path)

            logger.info(f"    - 同步 {patient_id} 的 FA_to_T1 和 Lesion Mask 的几何头文件信息...")
            patient_lesion_ants = ants.copy_image_info(patient_fa_ants, patient_lesion_ants)

            raw_mask_file = os.path.join(self.dirs['reg_to_con'], f"{patient_id}_lesion_raw.nii.gz")
            ants.image_write(patient_lesion_ants, raw_mask_file)

            for con_file in control_files:
                con_id_with_prefix = os.path.basename(con_file).split('.')[0]
                con_image_ants = ants.image_read(con_file)

                # 使用 FA_to_T1 进行配准
                registration = ants.registration(fixed=con_image_ants, moving=patient_fa_ants, type_of_transform="SyN")
                
                transformed_mask = ants.apply_transforms(fixed=con_image_ants, moving=patient_lesion_ants,
                                                         transformlist=registration['fwdtransforms'],
                                                         interpolator='nearestNeighbor')
                output_mask_file = os.path.join(self.dirs['reg_to_con'],
                                                f"{patient_id}_lesion_to_{con_id_with_prefix}.nii.gz")
                ants.image_write(transformed_mask, output_mask_file)

                transformed_fa = ants.apply_transforms(fixed=con_image_ants, moving=patient_fa_ants,
                                                       transformlist=registration['fwdtransforms'],
                                                       interpolator='linear')
                output_fa_file = os.path.join(self.dirs['reg_t1_to_con_qc'],
                                              f"{patient_id}_fa_to_{con_id_with_prefix}.nii.gz")
                ants.image_write(transformed_fa, output_fa_file)

                bedpostx_subject_id = con_id_with_prefix.split('_')[0]
                registered_lesion_info.append(
                    {'patient_id': patient_id, 'path': output_mask_file, 'base_id': bedpostx_subject_id})

        if not registered_lesion_info: raise PipelineStepEmptyError("Discon-Step 1 未生成配准文件。")
        return registered_lesion_info

    def _step2_probtrackx(self, registered_lesion_info):
        logger.info("Disconnection-Step 2: 运行纤维束追踪...")
        bedpostx_dir = self.config['disconnection_inputs']['bedpostx_dir']
        successful_probtrackx_outputs = 0

        for item in registered_lesion_info:
            lesion_mask_path, con_base_id = item['path'], item['base_id']
            mask_volume = get_fsl_volume(lesion_mask_path)
            if mask_volume == 0: continue

            subject_bedpostx_dir = os.path.join(bedpostx_dir, f"{con_base_id}.bedpostX")
            base_name = os.path.basename(lesion_mask_path).rsplit('.', 2)[0]
            output_dir = os.path.join(self.dirs['probtrackx'], f"{base_name}_output")
            os.makedirs(output_dir, exist_ok=True)

            cmd = (
                f"probtrackx2_gpu -x '{lesion_mask_path}' -l --onewaycondition -c 0.15 -S 2000 --steplength=0.15 -P 5000 "
                f"--fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd "
                f"-s '{os.path.join(subject_bedpostx_dir, 'merged')}' -m '{os.path.join(subject_bedpostx_dir, 'nodif_brain_mask')}' --dir='{output_dir}'")
            success, msg, _ = run_command_in_shell(cmd)
            if success: successful_probtrackx_outputs += 1

        if successful_probtrackx_outputs == 0: raise PipelineStepEmptyError("Discon-Step 2 未能生成 fdt_paths。")
        return self.dirs['probtrackx']

    def _intermediate_processing_chain(self, probtrackx_base_dir):
        logger.info("Disconnection-Step 3-5: 中间处理 (Python版)...")
        probtrackx_dirs = glob.glob(os.path.join(probtrackx_base_dir, '*_output'))
        binarized_fdt_paths = []
        for pdir in probtrackx_dirs:
            fdt_path_file = os.path.join(pdir, 'fdt_paths.nii.gz')
            if os.path.exists(fdt_path_file):
                base_name = os.path.basename(pdir).replace('_output', '')
                patient_id = base_name.split('_lesion_to_')[0]
                output_bin_file = os.path.join(self.dirs['fdt_paths'], f"{base_name}_fdt_bin.nii.gz")
                cmd = f"fslmaths '{fdt_path_file}' -bin '{output_bin_file}'"
                success, msg, _ = run_command_in_shell(cmd)
                if success: binarized_fdt_paths.append({'patient_id': patient_id, 'path': output_bin_file})

        registered_to_meanfa_files = []
        mat_dir, mean_fa_file = self.config['disconnection_inputs']['mat_dir'], self.config['disconnection_inputs'][
            'mean_fa_file']
        reference_image = sitk.ReadImage(mean_fa_file, sitk.sitkFloat32)

        for item in binarized_fdt_paths:
            fdt_bin_file, base_name = item['path'], os.path.basename(item['path']).replace('_fdt_bin.nii.gz', '')
            con_id_base = base_name.split('_to_')[-1].replace('registered_', '')
            ants_mat_file = os.path.join(mat_dir, f"{con_id_base}_fwd_1.mat")

            moving_image = sitk.ReadImage(fdt_bin_file, sitk.sitkUInt8)
            transform = sitk.ReadTransform(ants_mat_file)
            resampler = sitk.ResampleImageFilter()
            resampler.SetReferenceImage(reference_image)
            resampler.SetTransform(transform)
            resampler.SetInterpolator(sitk.sitkNearestNeighbor)
            resampler.SetDefaultPixelValue(0)
            resampler.SetOutputPixelType(sitk.sitkUInt8)
            resampled_image = resampler.Execute(moving_image)
            final_reg_file = os.path.join(self.dirs['reg_to_meanfa'], f"registered_{base_name}.nii.gz")
            sitk.WriteImage(resampled_image, final_reg_file)
            registered_to_meanfa_files.append({'patient_id': item['patient_id'], 'path': final_reg_file})

        grouped_files, patient_intersections = defaultdict(list), {}
        for item in registered_to_meanfa_files: grouped_files[item['patient_id']].append(item['path'])

        for patient_id, files in grouped_files.items():
            if not files: continue
            intersection_file = os.path.join(self.dirs['intersection'], f"{patient_id}_intersection.nii.gz")
            shutil.copy(files[0], intersection_file)
            for i in range(1, len(files)):
                run_command_in_shell(f"fslmaths '{intersection_file}' -mul '{files[i]}' '{intersection_file}'")
            if get_fsl_volume(intersection_file) > 0: patient_intersections[patient_id] = intersection_file

        return patient_intersections

    def _step6_calculate_results(self, patient_intersections, all_patient_ids):
        logger.info("Disconnection-Step 6: 计算最终结果...")
        all_results = []
        final_thr_map = self.config['disconnection_inputs']['final_thr_map']
        threshold_volume = get_fsl_volume(final_thr_map)

        for patient_id in all_patient_ids:
            if patient_id in patient_intersections:
                intersection_file = patient_intersections[patient_id]
                final_overlap_file = os.path.join(self.dirs['final_output'], f"{patient_id}_final_overlap.nii.gz")
                run_command_in_shell(f"fslmaths '{intersection_file}' -mul '{final_thr_map}' '{final_overlap_file}'")
                multiplied_volume = get_fsl_volume(final_overlap_file)
                percentage = (multiplied_volume / threshold_volume) * 100 if threshold_volume > 0 else 0
                all_results.append({
                    "Patient ID": patient_id, "Discon. Vol (vox)": multiplied_volume,
                    "Discon. %": f"{percentage:.2f}", "final_overlap_file": final_overlap_file
                })
            else:
                all_results.append({
                    "Patient ID": patient_id, "Discon. Vol (vox)": 0,
                    "Discon. %": "0.00 (No Intersection)", "final_overlap_file": "N/A"
                })
        return pd.DataFrame(all_results)
