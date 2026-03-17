import os
import glob
import shutil
import pandas as pd
import ants
import SimpleITK as sitk
from collections import defaultdict
from utils import run_command_in_shell, get_fsl_volume, PipelineStepEmptyError

class DisconnectionPipeline:
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.dirs = {}

    def prepare_dirs(self):
        work_dir = self.config['work_dir']
        self.dirs['reg_to_con'] = os.path.join(work_dir, 'discon_1_masks_to_controls')
        self.dirs['probtrackx'] = os.path.join(work_dir, 'discon_2_probtrackx')
        self.dirs['fdt_paths'] = os.path.join(work_dir, 'discon_3_fdt_bin')
        self.dirs['reg_to_meanfa'] = os.path.join(work_dir, 'discon_4_to_meanFA')
        self.dirs['intersection'] = os.path.join(work_dir, 'discon_5_intersections')
        self.dirs['final_output'] = os.path.join(work_dir, 'discon_6_final_outputs')
        for d in self.dirs.values(): os.makedirs(d, exist_ok=True)

    def run(self):
        self.logger.info("开始 Disconnection 分析...")
        self.prepare_dirs()
        reg_lesions = self.step1_registration()
        prob_dir = self.step2_probtrackx(reg_lesions)
        intersections = self.step3_5_intermediate(prob_dir)
        patient_ids = [p['id'] for p in self.config['patients']]
        df = self.step6_calculate(intersections, patient_ids)
        return df

    def step1_registration(self):
        self.logger.info("Discon-Step 1: FA to Controls FA Registration...")
        registered_info = []
        controls_dir = self.config['controls_FA_to_T1_dir']
        control_files = sorted(glob.glob(os.path.join(controls_dir, '*.nii.gz')))

        for p_info in self.config['patients']:
            p_id, p_fa, p_lesion = p_info['id'], p_info['fa_t1'], p_info['lesion']
            pat_fa_ants = ants.image_read(p_fa)
            pat_lesion_ants = ants.image_read(p_lesion)
            pat_lesion_ants = ants.copy_image_info(pat_fa_ants, pat_lesion_ants)

            for con_file in control_files:
                con_id = os.path.basename(con_file).split('.')[0]
                con_fa_ants = ants.image_read(con_file)

                # 使用 FA 配准 FA
                reg = ants.registration(fixed=con_fa_ants, moving=pat_fa_ants, type_of_transform="SyN")
                transformed_mask = ants.apply_transforms(
                    fixed=con_fa_ants, moving=pat_lesion_ants,
                    transformlist=reg['fwdtransforms'], interpolator='nearestNeighbor'
                )
                
                out_mask = os.path.join(self.dirs['reg_to_con'], f"{p_id}_lesion_to_{con_id}.nii.gz")
                ants.image_write(transformed_mask, out_mask)
                bedpostx_id = con_id.split('_')[0]
                registered_info.append({'patient_id': p_id, 'path': out_mask, 'base_id': bedpostx_id})

        return registered_info

    def step2_probtrackx(self, registered_info):
        self.logger.info("Discon-Step 2: Probtrackx2...")
        bedpostx_dir = self.config['bedpostx_dir']
        for item in registered_info:
            mask_path, base_id = item['path'], item['base_id']
            if get_fsl_volume(mask_path, self.logger) == 0: continue
            
            sub_bed_dir = os.path.join(bedpostx_dir, f"{base_id}.bedpostX")
            out_dir = os.path.join(self.dirs['probtrackx'], f"{os.path.basename(mask_path).split('.')[0]}_output")
            os.makedirs(out_dir, exist_ok=True)

            cmd = (f"probtrackx2_gpu -x '{mask_path}' -l --onewaycondition -c 0.15 -S 2000 --steplength=0.15 -P 5000 "
                   f"--fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd "
                   f"-s '{os.path.join(sub_bed_dir, 'merged')}' -m '{os.path.join(sub_bed_dir, 'nodif_brain_mask')}' --dir='{out_dir}'")
            run_command_in_shell(cmd, self.logger)
        return self.dirs['probtrackx']

    def step3_5_intermediate(self, prob_dir):
        self.logger.info("Discon-Step 3-5: Binarize, Transform & Intersect...")
        # 二值化
        bin_files = []
        for pdir in glob.glob(os.path.join(prob_dir, '*_output')):
            fdt = os.path.join(pdir, 'fdt_paths.nii.gz')
            if os.path.exists(fdt):
                base = os.path.basename(pdir).replace('_output', '')
                p_id = base.split('_lesion_to_')[0]
                out_bin = os.path.join(self.dirs['fdt_paths'], f"{base}_bin.nii.gz")
                run_command_in_shell(f"fslmaths '{fdt}' -bin '{out_bin}'", self.logger)
                bin_files.append({'p_id': p_id, 'path': out_bin})

        # SimpleITK 变换到 mean_fa
        reg_files = []
        ref_img = sitk.ReadImage(self.config['mean_fa'], sitk.sitkFloat32)
        for item in bin_files:
            base = os.path.basename(item['path']).replace('_bin.nii.gz', '')
            con_id = base.split('_to_')[-1]
            mat_file = os.path.join(self.config['mat_dir'], f"{con_id}_fwd_1.mat")
            
            mov_img = sitk.ReadImage(item['path'], sitk.sitkUInt8)
            transform = sitk.ReadTransform(mat_file)
            resampler = sitk.ResampleImageFilter()
            resampler.SetReferenceImage(ref_img)
            resampler.SetTransform(transform)
            resampler.SetInterpolator(sitk.sitkNearestNeighbor)
            resampled = resampler.Execute(mov_img)
            
            out_reg = os.path.join(self.dirs['reg_to_meanfa'], f"{base}_meanfa.nii.gz")
            sitk.WriteImage(resampled, out_reg)
            reg_files.append({'p_id': item['p_id'], 'path': out_reg})

        # 取交集
        intersections = {}
        grouped = defaultdict(list)
        for it in reg_files: grouped[it['p_id']].append(it['path'])
        
        for p_id, files in grouped.items():
            if not files: continue
            out_intersect = os.path.join(self.dirs['intersection'], f"{p_id}_intersect.nii.gz")
            shutil.copy(files[0], out_intersect)
            for i in range(1, len(files)):
                run_command_in_shell(f"fslmaths '{out_intersect}' -mul '{files[i]}' '{out_intersect}'", self.logger)
            if get_fsl_volume(out_intersect, self.logger) > 0:
                intersections[p_id] = out_intersect
        return intersections

    def step6_calculate(self, intersections, all_pids):
        self.logger.info("Discon-Step 6: 计算结果...")
        thr_map = self.config['thr_map']
        thr_vol = get_fsl_volume(thr_map, self.logger)
        results = []

        for p_id in all_pids:
            if p_id in intersections:
                final_overlap = os.path.join(self.dirs['final_output'], f"{p_id}_overlap.nii.gz")
                run_command_in_shell(f"fslmaths '{intersections[p_id]}' -mul '{thr_map}' '{final_overlap}'", self.logger)
                vol = get_fsl_volume(final_overlap, self.logger)
                pct = (vol / thr_vol) * 100 if thr_vol > 0 else 0
                results.append({"Patient ID": p_id, "Discon. Vol (vox)": vol, "Discon. %": f"{pct:.2f}"})
            else:
                results.append({"Patient ID": p_id, "Discon. Vol (vox)": 0, "Discon. %": "0.00"})
        return pd.DataFrame(results)
