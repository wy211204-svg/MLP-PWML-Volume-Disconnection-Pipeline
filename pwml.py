import os
import glob
import shutil
import re
import pandas as pd
from utils import run_command_in_shell, get_fsl_volume, PipelineStepEmptyError

class PWMLPipeline:
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.dirs = {}

    def prepare_dirs(self):
        work_dir = self.config['work_dir']
        self.dirs['reg'] = os.path.join(work_dir, 'pwml_1_reg')
        self.dirs['warp'] = os.path.join(work_dir, 'pwml_2_warped_atlas')
        self.dirs['results'] = os.path.join(work_dir, 'pwml_3_results')
        self.dirs['cluster'] = os.path.join(work_dir, 'pwml_4_cluster')
        for d in self.dirs.values(): os.makedirs(d, exist_ok=True)

    def run(self):
        self.logger.info("开始 PWML 分析...")
        self.prepare_dirs()
        results = []
        for p_info in self.config['patients']:
            trans = self.step1_register(p_info)
            atlas = self.step2_warp_atlas(p_info, trans)
            overlap = self.step3_calculate_overlap(p_info, atlas['bin'])
            msf = self.step4_cluster(p_info, atlas['bin'])
            fibers = self.step5_parse(p_info['id'], msf)
            results.append({**overlap, **fibers})
        return pd.DataFrame(results)

    def step1_register(self, p_info):
        p_id, p_fa = p_info['id'], p_info['fa_t1']
        mean_fa = self.config['mean_fa']
        self.logger.info(f"PWML-Step 1: FA to Mean_FA ({p_id})")
        
        reg_dir = os.path.join(self.dirs['reg'], p_id)
        os.makedirs(reg_dir, exist_ok=True)
        base = os.path.join(reg_dir, p_id)
        l6_out, l12_out, nonl_out = f"{base}_6dof.nii.gz", f"{base}_12dof.nii.gz", f"{base}_nonl.nii.gz"
        l6_mat, l12_mat, nonl_warp = f"{base}_6dof.mat", f"{base}_12dof.mat", f"{base}_nonl_warp.nii.gz"
        
        run_command_in_shell(f"flirt -in '{p_fa}' -ref '{mean_fa}' -dof 6 -cost corratio -out '{l6_out}' -omat '{l6_mat}'", self.logger)
        run_command_in_shell(f"flirt -in '{l6_out}' -ref '{mean_fa}' -dof 12 -cost corratio -out '{l12_out}' -omat '{l12_mat}'", self.logger)
        run_command_in_shell(f"fnirt --in='{l12_out}' --ref='{mean_fa}' --cout='{nonl_warp}' --iout='{nonl_out}' --config=FA_2_FMRIB58_1mm", self.logger)
        
        inv_l6, inv_l12, inv_warp = f"{base}_inv_6.mat", f"{base}_inv_12.mat", f"{base}_inv_warp.nii.gz"
        run_command_in_shell(f"convert_xfm -omat '{inv_l6}' -inverse '{l6_mat}'", self.logger)
        run_command_in_shell(f"convert_xfm -omat '{inv_l12}' -inverse '{l12_mat}'", self.logger)
        run_command_in_shell(f"invwarp --ref='{l12_out}' --warp='{nonl_warp}' --out='{inv_warp}'", self.logger)
        
        return {'inv_l6': inv_l6, 'inv_l12': inv_l12, 'inv_warp': inv_warp, 'ref_fa': p_fa, 'ref_l6': l6_out, 'ref_l12': l12_out}

    def step2_warp_atlas(self, p_info, trans):
        p_id = p_info['id']
        jhu = self.config['jhu_atlas']
        warp_dir = os.path.join(self.dirs['warp'], p_id)
        os.makedirs(warp_dir, exist_ok=True)
        tmp = os.path.join(warp_dir, "tmp")
        os.makedirs(tmp, exist_ok=True)
        
        run_command_in_shell(f"fslsplit '{jhu}' '{os.path.join(tmp, 'vol_')}' -t", self.logger)
        splits = sorted(glob.glob(os.path.join(tmp, 'vol_*.nii.gz')))
        warped = []
        for f in splits:
            base = os.path.basename(f).split('.')[0]
            inv1, inv2, inv3 = os.path.join(tmp, f'i1_{base}.nii.gz'), os.path.join(tmp, f'i2_{base}.nii.gz'), os.path.join(tmp, f'i3_{base}.nii.gz')
            run_command_in_shell(f"applywarp -i '{f}' -o '{inv1}' -r '{trans['ref_l12']}' -w '{trans['inv_warp']}' --interp=nn", self.logger)
            run_command_in_shell(f"flirt -in '{inv1}' -ref '{trans['ref_l6']}' -out '{inv2}' -applyxfm -init '{trans['inv_l12']}' -interp nearestneighbour", self.logger)
            run_command_in_shell(f"flirt -in '{inv2}' -ref '{trans['ref_fa']}' -out '{inv3}' -applyxfm -init '{trans['inv_l6']}' -interp nearestneighbour", self.logger)
            warped.append(inv3)
            
        prob = os.path.join(warp_dir, f'{p_id}_JHU_prob.nii.gz')
        run_command_in_shell(f"fslmerge -t '{prob}' {' '.join(warped)}", self.logger)
        bin_atlas = os.path.join(warp_dir, f'{p_id}_JHU_bin.nii.gz')
        run_command_in_shell(f"fslmaths '{prob}' -thr 15 -bin '{bin_atlas}'", self.logger)
        shutil.rmtree(tmp)
        return {'bin': bin_atlas}

    def step3_calculate_overlap(self, p_info, jhu_bin):
        p_id, p_lesion = p_info['id'], p_info['lesion']
        out = os.path.join(self.dirs['results'], f"{p_id}_lesion_in_wm.nii.gz")
        run_command_in_shell(f"fslmaths '{p_lesion}' -mul '{jhu_bin}' '{out}'", self.logger)
        l_vol = get_fsl_volume(out, self.logger)
        t_vol = get_fsl_volume(jhu_bin, self.logger)
        pct = (l_vol / t_vol * 100) if t_vol > 0 else 0
        return {"Patient ID": p_id, "PWML (%)": f"{pct:.2f}"}

    def step4_cluster(self, p_info, jhu_bin):
        p_id, p_lesion = p_info['id'], p_info['lesion']
        c_dir = os.path.join(self.dirs['cluster'], p_id)
        os.makedirs(c_dir, exist_ok=True)
        mul_out = os.path.join(c_dir, f"{p_id}_lesion_in_JHU.nii.gz")
        run_command_in_shell(f"fslmaths '{p_lesion}' -mul '{jhu_bin}' '{mul_out}'", self.logger)
        
        tmp = os.path.join(c_dir, "tmp")
        os.makedirs(tmp, exist_ok=True)
        run_command_in_shell(f"fslsplit '{mul_out}' '{os.path.join(tmp, 'vol_')}' -t", self.logger)
        msf = os.path.join(c_dir, f"{p_id}_clusters.msf")
        
        with open(msf, 'w') as f_msf:
            for f in sorted(glob.glob(os.path.join(tmp, 'vol_*.nii.gz'))):
                v_name = os.path.basename(f)
                f_msf.write(f"processing {v_name}\n")
                _, _, out = run_command_in_shell(f"cluster -i '{f}' -t 0.95 --mm", self.logger)
                f_msf.write(out + "\n")
        shutil.rmtree(tmp)
        return msf

    def step5_parse(self, p_id, msf):
        # 需求: PWML 的纤维束增加到4个 10,14,15,18，对应 left_IFOF, left_SLF, right_SLF, left_AF
        targets = {
            'left_IFOF_vol_mm3': "vol_0010",
            'left_SLF_vol_mm3': "vol_0014",
            'right_SLF_vol_mm3': "vol_0015",
            'left_AF_vol_mm3': "vol_0018"
        }
        vols = {k: 0 for k in targets}
        curr_k = None
        
        if os.path.exists(msf):
            with open(msf, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("processing"):
                        curr_k = next((k for k, v in targets.items() if v in line), None)
                        continue
                    if curr_k and re.match(r'^\d', line):
                        cols = line.split()
                        if len(cols) > 1: vols[curr_k] += int(cols[1])
                        
        conv = 0.8789
        return {k: f"{(v * conv):.2f}" for k, v in vols.items()}
