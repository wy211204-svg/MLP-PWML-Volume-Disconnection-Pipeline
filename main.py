import os
import argparse
import pandas as pd
from discon import DisconnectionPipeline
from pwml import PWMLPipeline
from utils import logger


def setup_directories(work_dir):
    dirs = {
        'work': work_dir,
        'reg_to_con': os.path.join(work_dir, 'discon_1_all_registered_masks'),
        'reg_to_con_qc': os.path.join(work_dir, 'discon_1a_registered_for_QC'),
        'probtrackx': os.path.join(work_dir, 'discon_2_all_probtrackx_outputs'),
        'fdt_paths': os.path.join(work_dir, 'discon_3_all_fdt_paths_binarized'),
        'reg_to_meanfa': os.path.join(work_dir, 'discon_4_all_registered_to_meanFA'),
        'intersection': os.path.join(work_dir, 'discon_5_all_intersections'),
        'final_output': os.path.join(work_dir, 'discon_6_final_outputs'),
        'pwml_reg': os.path.join(work_dir, 'pwml_1_registrations'),
        'pwml_warped_atlas': os.path.join(work_dir, 'pwml_2_warped_atlases'),
        'pwml_results': os.path.join(work_dir, 'pwml_3_results'),
        'pwml_cluster': os.path.join(work_dir, 'pwml_4_cluster_analysis')
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return dirs


def scan_patients(root_dir, t1_discon="FA_to_T1.nii.gz", t1_pwml="FA_to_T1.nii.gz", lesion="lesion_mask.nii.gz"):
    patients = []
    if not os.path.isdir(root_dir):
        logger.error(f"患者目录不存在: {root_dir}")
        return patients
    for pid in sorted(os.listdir(root_dir)):
        pdir = os.path.join(root_dir, pid)
        if os.path.isdir(pdir):
            pinfo = {
                'id': pid,
                't1_discon': os.path.join(pdir, t1_discon),
                't1_pwml': os.path.join(pdir, t1_pwml),
                'lesion': os.path.join(pdir, lesion)
            }
            if all(os.path.exists(v) for k, v in pinfo.items() if k != 'id'):
                patients.append(pinfo)
    return patients


def run_decision_tree(df):
    def apply_tree(row):
        try:
            discon = float(row.get('Discon_Score_Numeric', 0))
            left_ifof = float(row.get('Left_IFOF_vol_mm3', 0))
            right_slf = float(row.get('Right_SLF_vol_mm3', 0))
        except (ValueError, TypeError):
            return "Error"
        if discon < 0.63: return "MDI_normal"
        if left_ifof >= 0.44: return "MDI_delay"
        if right_slf >= 11.05: return "MDI_delay"
        return "MDI_normal"

    df['Decision_Tree_Prediction'] = df.apply(apply_tree, axis=1)
    return df


def main():
    parser = argparse.ArgumentParser(description="Disconnection & PWML Pipeline CLI")
    parser.add_argument("--patients_dir", required=True, help="患者数据的根目录")
    parser.add_argument("--work_dir", required=True, help="输出的工作目录")
    parser.add_argument("--controls_FA_to_T1_dir", required=True, help="对照组TFA配准到T1图像目录")
    parser.add_argument("--bedpostx_dir", required=True, help="BedpostX 数据目录")
    parser.add_argument("--mean_fa", required=True, help="Mean FA 文件路径")
    parser.add_argument("--mat_dir", required=True, help="转换矩阵文件目录")
    parser.add_argument("--thr_map", required=True, help="二值化失连接目标图文件路径")
    parser.add_argument("--jhu_atlas", required=True, help="JHU 概率图谱文件路径")

    parser.add_argument("--skip_discon", action="store_true", help="跳过 Disconnection 流程")
    parser.add_argument("--skip_pwml", action="store_true", help="跳过 PWML 流程")

    args = parser.parse_args()

    patients = scan_patients(args.patients_dir)
    if not patients:
        logger.error("未找到符合要求的患者数据，请检查目录和文件名！")
        return

    config = {
        'patients_to_process': patients,
        'disconnection_inputs': {
            'controls_t1_dir': args.controls_t1_dir,
            'bedpostx_dir': args.bedpostx_dir,
            'mean_fa_file': args.mean_fa,
            'mat_dir': args.mat_dir,
            'final_thr_map': args.thr_map,
        },
        'pwml_inputs': {
            'jhu_atlas_file': args.jhu_atlas
        }
    }

    dirs = setup_directories(args.work_dir)
    discon_df, pwml_df = None, None

    if not args.skip_discon:
        discon_df = DisconnectionPipeline(config, dirs).run()
    if not args.skip_pwml:
        pwml_df = PWMLPipeline(config, dirs).run()

    # 合并结果
    final_df = pd.DataFrame([p['id'] for p in patients], columns=["Patient ID"])
    if discon_df is not None and not discon_df.empty:
        cols_to_use = [c for c in discon_df.columns if c in ["Patient ID", "Discon. %", "Discon. Vol (vox)"]]
        final_df = pd.merge(final_df, discon_df[cols_to_use], on="Patient ID", how="left")
        
    if pwml_df is not None and not pwml_df.empty:
        final_df = pd.merge(final_df, pwml_df, on="Patient ID", how="left")


    out_csv = os.path.join(dirs['work'], 'final_results.csv')
    final_df.to_csv(out_csv, index=False)
    logger.info(f"所有流程执行完毕！最终结果保存在: {out_csv}")


if __name__ == "__main__":
    main()
