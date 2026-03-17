import os
import argparse
import pandas as pd
from utils import setup_logger
from discon_pipeline import DisconnectionPipeline
from pwml_pipeline import PWMLPipeline

def scan_patients(patients_dir, logger):
    patients = []
    if not os.path.exists(patients_dir):
        logger.error(f"找不到患者目录: {patients_dir}")
        return patients
        
    for p_id in sorted(os.listdir(patients_dir)):
        p_path = os.path.join(patients_dir, p_id)
        if os.path.isdir(p_path):
            fa_t1 = os.path.join(p_path, "FA_to_T1.nii.gz")
            lesion = os.path.join(p_path, "lesion_mask.nii.gz")
            if os.path.exists(fa_t1) and os.path.exists(lesion):
                patients.append({'id': p_id, 'fa_t1': fa_t1, 'lesion': lesion})
            else:
                logger.warning(f"跳过 {p_id}：未找到 FA_to_T1.nii.gz 或 lesion_mask.nii.gz")
    
    logger.info(f"扫描完成，找到 {len(patients)} 位患者数据。")
    return patients

def main():
    parser = argparse.ArgumentParser(description="Disconnection & PWML CLI Pipeline")
    parser.add_argument('--patients_dir', required=True, help="Root directory of patient data")
    parser.add_argument('--work_dir', required=True, help="Output workspace directory")
    parser.add_argument('--controls_FA_to_T1_dir', required=True, help="Controls FA to T1 directory")
    parser.add_argument('--bedpostx_dir', required=True, help="BedpostX results directory")
    parser.add_argument('--mean_fa', required=True, help="Study-specific FA template")
    parser.add_argument('--mat_dir', required=True, help="ANTs transformation matrices directory")
    parser.add_argument('--thr_map', required=True, help="Binary disconnection target map")
    parser.add_argument('--jhu_atlas', required=True, help="JHU white matter atlas")
    
    args = parser.parse_args()
    logger = setup_logger(args.work_dir)
    
    # 扫描患者
    patients = scan_patients(args.patients_dir, logger)
    if not patients:
        logger.error("未找到有效的患者数据，程序退出。")
        return

    config = vars(args)
    config['patients'] = patients

    df_final = pd.DataFrame([p['id'] for p in patients], columns=["Patient ID"])

    # 1. 运行 Disconnection Pipeline
    try:
        discon = DisconnectionPipeline(config, logger)
        df_discon = discon.run()
        df_final = pd.merge(df_final, df_discon, on="Patient ID", how="left")
    except Exception as e:
        logger.error(f"Disconnection 流程出错: {e}", exc_info=True)

    # 2. 运行 PWML Pipeline
    try:
        pwml = PWMLPipeline(config, logger)
        df_pwml = pwml.run()
        df_final = pd.merge(df_final, df_pwml, on="Patient ID", how="left")
    except Exception as e:
        logger.error(f"PWML 流程出错: {e}", exc_info=True)

    # 保存最终结果
    out_csv = os.path.join(args.work_dir, "final_results.csv")
    df_final.to_csv(out_csv, index=False)
    logger.info(f"🎉 全部流程执行完毕！最终结果已保存至: {out_csv}")
    print(df_final.to_string())

if __name__ == "__main__":
    main()
