import os
import argparse
import pandas as pd
from utils import setup_logger
from discon_pipeline import DisconnectionPipeline
from pwml_pipeline import PWMLPipeline

def scan_patients(patients_dir, logger):
    patients = []
    if not os.path.exists(patients_dir):
        logger.error(f"Patient directory not found: {patients_dir}")
        return patients
        
    for p_id in sorted(os.listdir(patients_dir)):
        p_path = os.path.join(patients_dir, p_id)
        if os.path.isdir(p_path):
            # Updated file name to T1.nii.gz
            t1 = os.path.join(p_path, "T1.nii.gz")
            lesion = os.path.join(p_path, "lesion_mask.nii.gz")
            if os.path.exists(t1) and os.path.exists(lesion):
                # Updated key name from 'fa_t1' to 't1'
                patients.append({'id': p_id, 't1': t1, 'lesion': lesion})
            else:
                logger.warning(f"Skipping {p_id}: T1.nii.gz or lesion_mask.nii.gz not found")
    
    logger.info(f"Scan complete. Found {len(patients)} patients with valid data.")
    return patients

def main():
    parser = argparse.ArgumentParser(description="Disconnection & PWML CLI Pipeline")
    parser.add_argument('--patients_dir', required=True, help="Root directory of patient data")
    parser.add_argument('--work_dir', required=True, help="Output workspace directory")
    # Updated to match --controls_T1_to_FA_dir
    parser.add_argument('--controls_T1_to_FA_dir', required=True, help="Directory containing control group T1 maps which registered to FA images")
    parser.add_argument('--bedpostx_dir', required=True, help="BedpostX results directory")
    # Updated to match --JHU_T1
    parser.add_argument('--JHU_T1', required=True, help="Study-specific FA template from JHU T1 template")
    parser.add_argument('--mat_dir', required=True, help="Directory containing ANTs transformation matrices")
    # Updated to match --JHU_atlas
    parser.add_argument('--JHU_atlas', required=True, help="JHU white matter atlas")
    
    args = parser.parse_args()
    logger = setup_logger(args.work_dir)
    
    # Scan for patients
    patients = scan_patients(args.patients_dir, logger)
    if not patients:
        logger.error("No valid patient data found. Exiting program.")
        return

    config = vars(args)
    config['patients'] = patients

    df_final = pd.DataFrame([p['id'] for p in patients], columns=["Patient ID"])

    # 1. Run Disconnection Pipeline
    try:
        discon = DisconnectionPipeline(config, logger)
        df_discon = discon.run()
        df_final = pd.merge(df_final, df_discon, on="Patient ID", how="left")
    except Exception as e:
        logger.error(f"Error in Disconnection pipeline: {e}", exc_info=True)

    # 2. Run PWML Pipeline
    try:
        pwml = PWMLPipeline(config, logger)
        df_pwml = pwml.run()
        df_final = pd.merge(df_final, df_pwml, on="Patient ID", how="left")
    except Exception as e:
        logger.error(f"Error in PWML pipeline: {e}", exc_info=True)

    # Save final results
    out_csv = os.path.join(args.work_dir, "final_results.csv")
    df_final.to_csv(out_csv, index=False)
    logger.info(f"All pipelines completed successfully! Final results saved to: {out_csv}")
    print(df_final.to_string())

if __name__ == "__main__":
    main()
