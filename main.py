import argparse
import os

from discon import DisconnectionPipeline
from utils import setup_logger


def scan_patients(patients_dir, logger):
    """Find patients containing both T1.nii.gz and lesion_mask.nii.gz."""

    patients = []

    if not os.path.exists(patients_dir):
        logger.error(f"Patient directory not found: {patients_dir}")
        return patients

    for patient_id in sorted(os.listdir(patients_dir)):
        patient_dir = os.path.join(patients_dir, patient_id)
        if not os.path.isdir(patient_dir):
            continue

        t1_path = os.path.join(patient_dir, "T1.nii.gz")
        lesion_path = os.path.join(patient_dir, "lesion_mask.nii.gz")

        if os.path.exists(t1_path) and os.path.exists(lesion_path):
            patients.append(
                {
                    "id": patient_id,
                    "t1": t1_path,
                    "lesion": lesion_path,
                }
            )
        else:
            logger.warning(
                f"Skipping {patient_id}: T1.nii.gz or lesion_mask.nii.gz not found"
            )

    logger.info(
        f"Patient scan complete. Found {len(patients)} patients with valid data."
    )
    return patients


def build_parser():
    parser = argparse.ArgumentParser(
        description="Structural Disconnectome Mapping Pipeline"
    )

    parser.add_argument(
        "--patients_dir",
        required=True,
        help="Root directory of patient data",
    )
    parser.add_argument(
        "--work_dir",
        required=True,
        help="Output workspace directory",
    )
    parser.add_argument(
        "--controls_T1_to_FA_dir",
        required=True,
        help=(
            "Directory containing control T1 images already aligned to the "
            "corresponding control FA/diffusion space"
        ),
    )
    parser.add_argument(
        "--bedpostx_dir",
        required=True,
        help="Directory containing <control_id>.bedpostX folders",
    )
    parser.add_argument(
        "--JHU_T1",
        required=True,
        help="Common JHU_T1 reference/template used for disconnectome maps",
    )
    parser.add_argument(
        "--mat_dir",
        required=True,
        help=(
            "Directory containing control-to-JHU SimpleITK-readable transforms "
            "named <control_id>_fwd_1.mat"
        ),
    )

    # Step-6 clinical grouping information is required to calculate the final Discon Score.
    parser.add_argument(
        "--clinical_groups_csv",
        required=True,
        help=(
            "CSV containing real clinical group assignments required for "
            "Step 6 Discon Score calculation"
        ),
    )
    parser.add_argument(
        "--patient_id_column",
        default="patient_id",
        help="Patient ID column in --clinical_groups_csv (default: patient_id)",
    )
    parser.add_argument(
        "--group_column",
        default="group",
        help="Group column in --clinical_groups_csv (default: group)",
    )
    parser.add_argument(
        "--unimpaired_label",
        default="unimpaired",
        help="Label representing the unimpaired group (default: unimpaired)",
    )
    parser.add_argument(
        "--delay_label",
        default="delay",
        help="Label representing the delay group (default: delay)",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    logger = setup_logger(args.work_dir)

    patients = scan_patients(args.patients_dir, logger)
    if not patients:
        logger.error("No valid patient data found. Exiting program.")
        return

    config = vars(args)
    config["patients"] = patients

    try:
        pipeline = DisconnectionPipeline(config, logger)
        df_final = pipeline.run()
    except Exception as exc:
        logger.error(
            f"Error in Structural Disconnection pipeline: {exc}",
            exc_info=True,
        )
        raise

    output_csv = os.path.join(args.work_dir, "final_results.csv")
    df_final.to_csv(output_csv, index=False)

    logger.info(
        "Structural Disconnection pipeline completed successfully. "
        f"Final results saved to: {output_csv}"
    )
    print(df_final.to_string(index=False))


if __name__ == "__main__":
    main()
