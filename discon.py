import os
import glob
import pandas as pd
import numpy as np
import nibabel as nib
import ants
import SimpleITK as sitk
from collections import defaultdict
from utils import run_command_in_shell, get_fsl_volume

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
        self.dirs['patient_avg'] = os.path.join(work_dir, 'discon_5_patient_averages')
        self.dirs['final_output'] = os.path.join(work_dir, 'discon_6_final_outputs')
        for d in self.dirs.values(): os.makedirs(d, exist_ok=True)

    def run(self):
        self.logger.info("Starting Disconnection analysis...")
        self.prepare_dirs()
        
        reg_lesions = self.step1_registration()
        prob_dir = self.step2_probtrackx(reg_lesions)
        patient_avg_maps = self.step3_5_stats_and_average(prob_dir)
        
        patient_ids = [p['id'] for p in self.config['patients']]
        df = self.step6_group_calculate(patient_avg_maps, patient_ids)
        return df

    def step1_registration(self):
        self.logger.info("Discon-Step 1: FA to Controls FA Registration...")
        registered_info = []
        # Name updated: controls_FA_to_T1_dir -> controls_T1_to_FA_dir
        controls_dir = self.config['controls_T1_to_FA_dir']
        control_files = sorted(glob.glob(os.path.join(controls_dir, '*.nii.gz')))

        for p_info in self.config['patients']:
            p_id, p_fa, p_lesion = p_info['id'], p_info['fa_t1'], p_info['lesion']
            pat_fa_ants = ants.image_read(p_fa)
            pat_lesion_ants = ants.image_read(p_lesion)
            pat_lesion_ants = ants.copy_image_info(pat_fa_ants, pat_lesion_ants)

            for con_file in control_files:
                con_id = os.path.basename(con_file).split('.')[0]
                con_fa_ants = ants.image_read(con_file)

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

    def step3_5_stats_and_average(self, prob_dir):
        self.logger.info("Discon-Step 3-5: Stats Thresholding, Transform & Patient Average...")
        
        # --- Step 3: Threshold based on mean + 2*std ---
        bin_files = []
        for pdir in glob.glob(os.path.join(prob_dir, '*_output')):
            fdt = os.path.join(pdir, 'fdt_paths.nii.gz')
            if os.path.exists(fdt):
                base = os.path.basename(pdir).replace('_output', '')
                p_id = base.split('_lesion_to_')[0]
                out_bin = os.path.join(self.dirs['fdt_paths'], f"{base}_bin.nii.gz")
                
                img = nib.load(fdt)
                data = img.get_fdata()
                
                nonzero_data = data[data > 0]
                if len(nonzero_data) == 0:
                    continue
                
                mean_val = np.mean(nonzero_data)
                std_val = np.std(nonzero_data)
                threshold = mean_val + 2 * std_val
                
                bin_data = (data >= threshold).astype(np.uint8)
                bin_img = nib.Nifti1Image(bin_data, img.affine, img.header)
                nib.save(bin_img, out_bin)
                
                bin_files.append({'p_id': p_id, 'path': out_bin})

        # --- Step 4: SimpleITK transform to JHU_T1 space ---
        reg_files = []
        # Name updated: mean_fa -> JHU_T1
        ref_img = sitk.ReadImage(self.config['JHU_T1'], sitk.sitkFloat32)
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
            
            # Output name updated to reflect JHU_T1 space
            out_reg = os.path.join(self.dirs['reg_to_meanfa'], f"{base}_JHU_T1.nii.gz")
            sitk.WriteImage(resampled, out_reg)
            reg_files.append({'p_id': item['p_id'], 'path': out_reg})

        # --- Step 5: Average maps per patient ---
        patient_avg_maps = {}
        grouped = defaultdict(list)
        for it in reg_files: grouped[it['p_id']].append(it['path'])
        
        for p_id, files in grouped.items():
            if not files: continue
            
            sum_data = None
            ref_nii = None
            num_controls = len(files)
            
            for f in files:
                img = nib.load(f)
                data = img.get_fdata().astype(np.float32)
                if sum_data is None:
                    sum_data = data
                    ref_nii = img
                else:
                    sum_data += data
            
            avg_data = sum_data / num_controls
            avg_img = nib.Nifti1Image(avg_data, ref_nii.affine, ref_nii.header)
            
            out_avg = os.path.join(self.dirs['patient_avg'], f"{p_id}_avg.nii.gz")
            nib.save(avg_img, out_avg)
            patient_avg_maps[p_id] = out_avg
            
        return patient_avg_maps

    def step6_group_calculate(self, patient_avg_maps, all_pids):
        self.logger.info("Discon-Step 6: Group-wise calculation of final scores...")
        
        # Split patients into groups based on clinical scores
        half_idx = max(1, len(all_pids)//2)
        unimpaired_pids = all_pids[:half_idx]
        delay_pids = all_pids[half_idx:]
        
        def calculate_group_average(pid_list):
            sum_data, count, ref_img = None, 0, None
            for p_id in pid_list:
                if p_id in patient_avg_maps:
                    img = nib.load(patient_avg_maps[p_id])
                    if sum_data is None:
                        sum_data = img.get_fdata().astype(np.float32)
                        ref_img = img
                    else:
                        sum_data += img.get_fdata()
                    count += 1
            if count == 0: return None, None
            return sum_data / count, ref_img

        unimpaired_avg_data, ref_nii = calculate_group_average(unimpaired_pids)
        delay_avg_data, _ = calculate_group_average(delay_pids)

        if unimpaired_avg_data is None or delay_avg_data is None:
            self.logger.error("Insufficient group data to generate threshold; setting scores to 0.")
            return pd.DataFrame([{"Patient ID": p, "Pat Bin Vol (vox)": 0, "Target Vol (vox)": 0, "Discon Score": "0.0000"} for p in all_pids])

        max_unimpaired_val = np.max(unimpaired_avg_data)
        self.logger.info(f"Using max value from unimpaired group as threshold: {max_unimpaired_val:.6f}")

        delay_binarized_data = (delay_avg_data >= max_unimpaired_val).astype(np.uint8)
        target_map_path = os.path.join(self.dirs['final_output'], "Target_DelayGroup_Binarized.nii.gz")
        nib.save(nib.Nifti1Image(delay_binarized_data, ref_nii.affine, ref_nii.header), target_map_path)
        
        target_volume = np.sum(delay_binarized_data > 0)

        results = []
        for p_id in all_pids:
            if p_id in patient_avg_maps and target_volume > 0:
                pat_avg_data = nib.load(patient_avg_maps[p_id]).get_fdata()
                pat_binarized = (pat_avg_data > 0).astype(np.uint8)
                pat_volume = np.sum(pat_binarized > 0)
                score = pat_volume / target_volume
                
                results.append({
                    "Patient ID": p_id, 
                    "Pat Bin Vol (vox)": pat_volume,
                    "Target Vol (vox)": target_volume,
                    "Discon Score": f"{score:.4f}"
                })
            else:
                results.append({
                    "Patient ID": p_id, 
                    "Pat Bin Vol (vox)": 0,
                    "Target Vol (vox)": target_volume,
                    "Discon Score": "0.0000"
                })

        return pd.DataFrame(results)
