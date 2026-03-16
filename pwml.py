import os
import glob
import shutil
import re
import pandas as pd
from utils import run_command_in_shell, get_fsl_volume, PipelineStepEmptyError, logger


class PWMLPipeline:
    def __init__(self, config, dirs):
        self.config = config
        self.dirs = dirs

    def run(self):
        logger.info("PWML: 开始分析...")
        all_results = []
        for p_info in self.config['patients_to_process']:
            logger.info(f"--- PWML: 开始处理患者 {p_info['id']} ---")
            transform_paths = self._step1_register_t1_to_fa(p_info)
            native_space_atlas = self._step2_warp_atlas_to_t1(p_info, transform_paths)
            pwml_result = self._step3_calculate_overlap(p_info, native_space_atlas['bin'])
            msf_file_path = self._step4_run_cluster_analysis(p_info, native_space_atlas, self.dirs['pwml_cluster'])
            decision_tree_results = self._step5_parse_cluster_results(p_info['id'], msf_file_path)
            combined_result = {**pwml_result, **decision_tree_results}
            all_results.append(combined_result)

        if not all_results: raise PipelineStepEmptyError("PWML pipeline did not produce any results.")
        return pd.DataFrame(all_results)

    def _step1_register_t1_to_fa(self, p_info):
        patient_id, patient_t1 = p_info['id'], p_info['t1_pwml']
        target_fa = self.config['pwml_inputs']['target_fa_file']
        fnirt_config_name = "FA_2_FMRIB58_1mm"

        reg_dir = os.path.join(self.dirs['pwml_reg'], patient_id)
        os.makedirs(reg_dir, exist_ok=True)
        base = os.path.join(reg_dir, patient_id)

        l6_out, l12_out, nonl_out = f"{base}_6dof.nii.gz", f"{base}_12dof.nii.gz", f"{base}_nonl.nii.gz"
        l6_mat, l12_mat, nonl_warp = f"{base}_6dof.mat", f"{base}_12dof.mat", f"{base}_nonl_warp.nii.gz"

        cmds = [
            f"flirt -in '{patient_t1}' -ref '{target_fa}' -dof 6 -cost corratio -out '{l6_out}' -omat '{l6_mat}'",
            f"flirt -in '{l6_out}' -ref '{target_fa}' -dof 12 -cost corratio -out '{l12_out}' -omat '{l12_mat}'",
            f"fnirt --in='{l12_out}' --ref='{target_fa}' --cout='{nonl_warp}' --iout='{nonl_out}' --config={fnirt_config_name}"
        ]
        for cmd in cmds: run_command_in_shell(cmd)

        inv_l6_mat, inv_l12_mat, inv_nonl_warp = f"{base}_inv_6dof.mat", f"{base}_inv_12dof.mat", f"{base}_inv_nonl_warp.nii.gz"
        cmds_inv = [
            f"convert_xfm -omat '{inv_l6_mat}' -inverse '{l6_mat}'",
            f"convert_xfm -omat '{inv_l12_mat}' -inverse '{l12_mat}'",
            f"invwarp --ref='{l12_out}' --warp='{nonl_warp}' --out='{inv_nonl_warp}'"
        ]
        for cmd in cmds_inv: run_command_in_shell(cmd)

        return {'inv_l6': inv_l6_mat, 'inv_l12': inv_l12_mat, 'inv_warp': inv_nonl_warp, 'ref_t1': patient_t1,
                'ref_l6': l6_out, 'ref_l12': l12_out}

    def _step2_warp_atlas_to_t1(self, p_info, transforms):
        patient_id = p_info['id']
        jhu_atlas = self.config['pwml_inputs']['jhu_atlas_file']
        warp_dir = os.path.join(self.dirs['pwml_warped_atlas'], patient_id)
        temp_split_dir = os.path.join(warp_dir, "temp_split_atlas")
        os.makedirs(temp_split_dir, exist_ok=True)

        run_command_in_shell(f"fslsplit '{jhu_atlas}' '{os.path.join(temp_split_dir, 'vol_')}' -t")
        split_files = sorted(glob.glob(os.path.join(temp_split_dir, 'vol_*.nii.gz')))
        final_warped_files = []

        for f_in in split_files:
            base_name = os.path.basename(f_in).split('.')[0]
            inv1, inv2, inv3 = [os.path.join(temp_split_dir, f'inv{i}_{base_name}.nii.gz') for i in (1, 2, 3)]
            final_warped_files.append(inv3)
            run_command_in_shell(
                f"applywarp -i '{f_in}' -o '{inv1}' -r '{transforms['ref_l12']}' -w '{transforms['inv_warp']}' --interp=nn")
            run_command_in_shell(
                f"flirt -in '{inv1}' -ref '{transforms['ref_l6']}' -out '{inv2}' -applyxfm -init '{transforms['inv_l12']}' -interp nearestneighbour")
            run_command_in_shell(
                f"flirt -in '{inv2}' -ref '{transforms['ref_t1']}' -out '{inv3}' -applyxfm -init '{transforms['inv_l6']}' -interp nearestneighbour")

        final_warped_atlas = os.path.join(warp_dir, f'{patient_id}_JHU_in_T1_space_prob.nii.gz')
        run_command_in_shell(f"fslmerge -t '{final_warped_atlas}' {' '.join(final_warped_files)}")
        shutil.rmtree(temp_split_dir)

        final_binned_atlas = os.path.join(warp_dir, f'{patient_id}_JHU_in_T1_space_bin.nii.gz')
        run_command_in_shell(f"fslmaths '{final_warped_atlas}' -thr 15 -bin '{final_binned_atlas}'")
        return {'prob': final_warped_atlas, 'bin': final_binned_atlas}

    def _step3_calculate_overlap(self, p_info, patient_wm_atlas_bin_path):
        patient_id, lesion_mask = p_info['id'], p_info['lesion']
        result_dir = os.path.join(self.dirs['pwml_results'], patient_id)
        os.makedirs(result_dir, exist_ok=True)
        lesion_in_wm_path = os.path.join(result_dir, f"{patient_id}_lesion_in_wm.nii.gz")
        run_command_in_shell(f"fslmaths '{lesion_mask}' -mul '{patient_wm_atlas_bin_path}' '{lesion_in_wm_path}'")

        lesioned_wm, total_wm = get_fsl_volume(lesion_in_wm_path), get_fsl_volume(patient_wm_atlas_bin_path)
        percentage = (lesioned_wm / total_wm) * 100 if total_wm > 0 else 0
        return {"Patient ID": patient_id, "Lesioned WM Vol (vox)": lesioned_wm, "Total WM Vol (vox)": total_wm,
                "PWML (%)": f"{percentage:.2f}"}

    def _step4_run_cluster_analysis(self, p_info, native_space_atlas, output_dir):
        patient_id, lesion_mask = p_info['id'], p_info['lesion']
        fsl_cluster_cmd = "cluster"  # Assuming it's in system PATH
        cluster_dir = os.path.join(output_dir, patient_id)
        os.makedirs(cluster_dir, exist_ok=True)

        lesion_in_bin_atlas = os.path.join(cluster_dir, f"{patient_id}_lesion_in_JHU_bin.nii.gz")
        run_command_in_shell(f"fslmaths '{lesion_mask}' -mul '{native_space_atlas['bin']}' '{lesion_in_bin_atlas}'")

        temp_split = os.path.join(cluster_dir, "temp_split")
        os.makedirs(temp_split, exist_ok=True)
        run_command_in_shell(f"fslsplit '{lesion_in_bin_atlas}' '{os.path.join(temp_split, 'vol_')}' -t")

        msf_file_path = os.path.join(cluster_dir, f"{patient_id}_clusters.msf")
        if os.path.exists(msf_file_path): os.remove(msf_file_path)

        for f in sorted(glob.glob(os.path.join(temp_split, 'vol_*.nii.gz'))):
            vol_name = os.path.basename(f)
            with open(msf_file_path, 'a') as msf_file:
                msf_file.write(f"processing {vol_name}\n")
            success, _, output = run_command_in_shell(f"'{fsl_cluster_cmd}' -i '{f}' -t 0.95 --mm")
            if success:
                with open(msf_file_path, 'a') as msf_file: msf_file.write(output + "\n")

        shutil.rmtree(temp_split)
        return msf_file_path

    def _step5_parse_cluster_results(self, patient_id, msf_file_path):
        target_fibers = {'Right_IFOF_vol_mm3': "vol_0015", 'Left_SLF_vol_mm3': "vol_0010"}
        fiber_volumes = {key: 0 for key in target_fibers}
        current_key = None

        if os.path.exists(msf_file_path):
            with open(msf_file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("processing"):
                        current_key = next((k for k, p in target_fibers.items() if f" {p}" in line), None)
                        continue
                    if current_key and re.match(r'^\d', line):
                        cols = line.split()
                        if len(cols) > 1: fiber_volumes[current_key] += int(cols[1])

        factor = 0.8789
        return {k: f"{(v * factor):.2f}" for k, v in fiber_volumes.items()}