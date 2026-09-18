import os
import glob
from collections import defaultdict

import ants
import nibabel as nib
import numpy as np
import pandas as pd
import SimpleITK as sitk

from utils import get_fsl_volume, run_command_in_shell


class DisconnectionPipeline:
    """Structural disconnectome pipeline.

    Pipeline overview
    -----------------
    Step 1: Register each patient's T1/lesion mask to each control space.
    Step 2: Run probtrackx2_gpu from the transformed lesion mask in each control.
    Step 3: Transform raw fdt_paths maps to the common JHU_T1 space.
            Per project specification, nearest-neighbor interpolation is retained.
    Step 4: Threshold each JHU-space fdt_paths map using mean(nonzero)+2*SD(nonzero)
            and generate a binary disconnectome map.
    Step 5: Average all binary control-derived maps for each patient to create a
            patient-level disconnectome probability/consistency map (0-1).
    Step 6: Use real clinical grouping data to construct a delay-group target
            and calculate each patient's Discon Score as:

                overlap(patient_disconnectome_binary, target) / target_volume

            Step 6 is required to obtain the final Discon Score.
    """

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.dirs = {}

    def prepare_dirs(self):
        work_dir = self.config["work_dir"]

        self.dirs["reg_to_con"] = os.path.join(
            work_dir, "discon_1_masks_to_controls"
        )
        self.dirs["probtrackx"] = os.path.join(
            work_dir, "discon_2_probtrackx"
        )
        self.dirs["fdt_to_jhu"] = os.path.join(
            work_dir, "discon_3_fdt_to_JHU"
        )
        self.dirs["binary_jhu"] = os.path.join(
            work_dir, "discon_4_binary_JHU"
        )
        self.dirs["patient_avg"] = os.path.join(
            work_dir, "discon_5_patient_disconnectome"
        )
        self.dirs["final_output"] = os.path.join(
            work_dir, "discon_6_clinical_outputs"
        )

        for directory in self.dirs.values():
            os.makedirs(directory, exist_ok=True)

    def run(self):
        self.logger.info("Starting Structural Disconnection analysis...")
        self.prepare_dirs()

        registered_lesions = self.step1_registration()
        prob_dir = self.step2_probtrackx(registered_lesions)
        jhu_fdt_files = self.step3_transform_fdt_to_jhu(prob_dir)
        binary_files = self.step4_threshold_in_jhu(jhu_fdt_files)
        patient_avg_maps = self.step5_patient_average(binary_files)

        # Step 5 provides the patient-level disconnectome maps used by Step 6.
        all_patient_ids = [p["id"] for p in self.config["patients"]]
        result_rows = []
        for patient_id in all_patient_ids:
            result_rows.append(
                {
                    "Patient ID": patient_id,
                    "Disconnectome Map": patient_avg_maps.get(patient_id, ""),
                }
            )
        df_results = pd.DataFrame(result_rows)

        # Step 6 is required for the final Discon Score.
        df_step6 = self.step6_group_calculate(patient_avg_maps)
        if df_step6 is None or df_step6.empty:
            raise RuntimeError("Step 6 did not produce Discon Scores.")

        df_results = pd.merge(
            df_results, df_step6, on="Patient ID", how="left"
        )

        return df_results

    def step1_registration(self):
        self.logger.info(
            "Discon-Step 1: Register patient T1/lesion masks to each control space..."
        )

        registered_info = []
        controls_dir = self.config["controls_T1_to_FA_dir"]
        control_files = sorted(glob.glob(os.path.join(controls_dir, "*.nii.gz")))

        if not control_files:
            raise RuntimeError(
                f"No control NIfTI files (*.nii.gz) found in: {controls_dir}"
            )

        for patient_info in self.config["patients"]:
            patient_id = patient_info["id"]
            patient_t1 = patient_info["t1"]
            patient_lesion = patient_info["lesion"]

            patient_t1_ants = ants.image_read(patient_t1)
            patient_lesion_ants = ants.image_read(patient_lesion)

            # Keep lesion-mask geometry consistent with the patient T1.
            patient_lesion_ants = ants.copy_image_info(
                patient_t1_ants, patient_lesion_ants
            )

            for control_file in control_files:
                control_id = os.path.basename(control_file).split(".")[0]
                control_ants = ants.image_read(control_file)

                self.logger.info(
                    f"Step 1: {patient_id} -> control {control_id}"
                )

                registration = ants.registration(
                    fixed=control_ants,
                    moving=patient_t1_ants,
                    type_of_transform="SyN",
                )

                transformed_mask = ants.apply_transforms(
                    fixed=control_ants,
                    moving=patient_lesion_ants,
                    transformlist=registration["fwdtransforms"],
                    interpolator="nearestNeighbor",
                )

                output_mask = os.path.join(
                    self.dirs["reg_to_con"],
                    f"{patient_id}_lesion_to_{control_id}.nii.gz",
                )
                ants.image_write(transformed_mask, output_mask)

                # Keep the original mapping rule used by the repository:
                # control filename prefix before '_' corresponds to <ID>.bedpostX.
                bedpostx_id = control_id.split("_")[0]

                registered_info.append(
                    {
                        "patient_id": patient_id,
                        "control_id": control_id,
                        "base_id": bedpostx_id,
                        "path": output_mask,
                    }
                )

        if not registered_info:
            raise RuntimeError("Step 1 produced no registered lesion masks.")

        return registered_info

    def step2_probtrackx(self, registered_info):
        self.logger.info("Discon-Step 2: Run probtrackx2_gpu...")

        bedpostx_dir = self.config["bedpostx_dir"]
        success_count = 0

        for item in registered_info:
            mask_path = item["path"]
            base_id = item["base_id"]
            patient_id = item["patient_id"]
            control_id = item["control_id"]

            if get_fsl_volume(mask_path, self.logger) == 0:
                self.logger.warning(
                    f"Skipping empty transformed lesion mask: {mask_path}"
                )
                continue

            sub_bed_dir = os.path.join(bedpostx_dir, f"{base_id}.bedpostX")
            merged_path = os.path.join(sub_bed_dir, "merged")
            nodif_mask = os.path.join(sub_bed_dir, "nodif_brain_mask")

            if not os.path.exists(sub_bed_dir):
                self.logger.warning(
                    f"BedpostX directory not found for control {control_id}: "
                    f"{sub_bed_dir}"
                )
                continue

            out_dir = os.path.join(
                self.dirs["probtrackx"],
                f"{patient_id}_lesion_to_{control_id}_output",
            )
            os.makedirs(out_dir, exist_ok=True)

            cmd = (
                f"probtrackx2_gpu -x '{mask_path}' -l --onewaycondition "
                f"-c 0.15 -S 2000 --steplength=0.15 -P 5000 "
                f"--fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 "
                f"--forcedir --opd "
                f"-s '{merged_path}' "
                f"-m '{nodif_mask}' "
                f"--dir='{out_dir}'"
            )

            ok, status, output = run_command_in_shell(cmd, self.logger)
            if not ok:
                self.logger.error(
                    f"probtrackx2_gpu failed for {patient_id} / {control_id}. "
                    f"{status}\n{output}"
                )
                continue

            fdt_path = os.path.join(out_dir, "fdt_paths.nii.gz")
            if not os.path.exists(fdt_path):
                self.logger.warning(
                    f"probtrackx2_gpu completed but fdt_paths.nii.gz was not found: "
                    f"{out_dir}"
                )
                continue

            success_count += 1

        if success_count == 0:
            raise RuntimeError("Step 2 produced no valid fdt_paths.nii.gz files.")

        self.logger.info(
            f"Discon-Step 2 completed: {success_count} valid probtrackx outputs."
        )
        return self.dirs["probtrackx"]

    def step3_transform_fdt_to_jhu(self, prob_dir):
        """Transform raw fdt_paths maps to JHU_T1 space.

        IMPORTANT: By project specification, nearest-neighbor interpolation is
        intentionally retained here even though fdt_paths is a continuous map.
        """

        self.logger.info(
            "Discon-Step 3: Transform raw fdt_paths maps to JHU_T1 space "
            "using nearest-neighbor interpolation..."
        )

        ref_img = sitk.ReadImage(self.config["JHU_T1"], sitk.sitkFloat32)
        transformed_files = []

        for prob_output_dir in sorted(
            glob.glob(os.path.join(prob_dir, "*_output"))
        ):
            fdt_path = os.path.join(prob_output_dir, "fdt_paths.nii.gz")
            if not os.path.exists(fdt_path):
                continue

            base = os.path.basename(prob_output_dir).replace("_output", "")
            if "_lesion_to_" not in base:
                self.logger.warning(
                    f"Unexpected probtrackx output directory name; skipping: {base}"
                )
                continue

            patient_id, control_id = base.split("_lesion_to_", 1)
            mat_file = os.path.join(
                self.config["mat_dir"], f"{control_id}_fwd_1.mat"
            )

            if not os.path.exists(mat_file):
                self.logger.warning(
                    f"Transform file not found for control {control_id}: {mat_file}"
                )
                continue

            moving_img = sitk.ReadImage(fdt_path, sitk.sitkFloat32)
            transform = sitk.ReadTransform(mat_file)

            resampler = sitk.ResampleImageFilter()
            resampler.SetReferenceImage(ref_img)
            resampler.SetTransform(transform)
            resampler.SetInterpolator(sitk.sitkNearestNeighbor)
            resampler.SetDefaultPixelValue(0.0)

            resampled = resampler.Execute(moving_img)

            output_path = os.path.join(
                self.dirs["fdt_to_jhu"],
                f"{base}_JHU_T1_fdt_paths.nii.gz",
            )
            sitk.WriteImage(resampled, output_path)

            transformed_files.append(
                {
                    "patient_id": patient_id,
                    "control_id": control_id,
                    "path": output_path,
                }
            )

        if not transformed_files:
            raise RuntimeError("Step 3 produced no JHU-space fdt_paths maps.")

        self.logger.info(
            f"Discon-Step 3 completed: {len(transformed_files)} maps transformed."
        )
        return transformed_files

    def step4_threshold_in_jhu(self, transformed_files):
        """Threshold JHU-space fdt_paths maps using mean(nonzero) + 2*SD."""

        self.logger.info(
            "Discon-Step 4: Threshold JHU-space fdt_paths maps using "
            "mean(nonzero) + 2*SD(nonzero)..."
        )

        binary_files = []
        threshold_records = []

        for item in transformed_files:
            image = nib.load(item["path"])
            data = image.get_fdata(dtype=np.float32)
            nonzero_data = data[data > 0]

            if nonzero_data.size == 0:
                self.logger.warning(
                    f"Skipping empty JHU-space fdt map: {item['path']}"
                )
                continue

            mean_val = float(np.mean(nonzero_data))
            std_val = float(np.std(nonzero_data))
            threshold = mean_val + 2.0 * std_val

            binary_data = (data >= threshold).astype(np.uint8)
            binary_voxels = int(np.count_nonzero(binary_data))

            header = image.header.copy()
            header.set_data_dtype(np.uint8)
            binary_img = nib.Nifti1Image(binary_data, image.affine, header)

            output_path = os.path.join(
                self.dirs["binary_jhu"],
                f"{item['patient_id']}_lesion_to_{item['control_id']}_JHU_T1_bin.nii.gz",
            )
            nib.save(binary_img, output_path)

            binary_files.append(
                {
                    "patient_id": item["patient_id"],
                    "control_id": item["control_id"],
                    "path": output_path,
                }
            )

            threshold_records.append(
                {
                    "Patient ID": item["patient_id"],
                    "Control ID": item["control_id"],
                    "Nonzero Voxels": int(nonzero_data.size),
                    "Mean Nonzero FDT": mean_val,
                    "SD Nonzero FDT": std_val,
                    "Threshold (Mean+2SD)": threshold,
                    "Binary Voxels": binary_voxels,
                }
            )

        if not binary_files:
            raise RuntimeError("Step 4 produced no binary disconnectome maps.")

        threshold_csv = os.path.join(
            self.dirs["binary_jhu"], "threshold_statistics.csv"
        )
        pd.DataFrame(threshold_records).to_csv(threshold_csv, index=False)
        self.logger.info(f"Threshold statistics saved to: {threshold_csv}")

        return binary_files

    def step5_patient_average(self, binary_files):
        """Average control-derived binary maps for each patient."""

        self.logger.info(
            "Discon-Step 5: Average binary maps across controls for each patient..."
        )

        grouped = defaultdict(list)
        for item in binary_files:
            grouped[item["patient_id"]].append(item["path"])

        patient_avg_maps = {}
        summary_rows = []

        for patient_id, files in grouped.items():
            if not files:
                continue

            sum_data = None
            ref_nii = None

            for file_path in files:
                image = nib.load(file_path)
                data = image.get_fdata(dtype=np.float32)

                if sum_data is None:
                    sum_data = data.copy()
                    ref_nii = image
                else:
                    if data.shape != sum_data.shape:
                        raise RuntimeError(
                            f"Shape mismatch in Step 5 for patient {patient_id}: "
                            f"{file_path}"
                        )
                    sum_data += data

            num_controls = len(files)
            avg_data = (sum_data / float(num_controls)).astype(np.float32)

            header = ref_nii.header.copy()
            header.set_data_dtype(np.float32)
            avg_img = nib.Nifti1Image(avg_data, ref_nii.affine, header)

            output_path = os.path.join(
                self.dirs["patient_avg"],
                f"{patient_id}_disconnectome_prob.nii.gz",
            )
            nib.save(avg_img, output_path)

            patient_avg_maps[patient_id] = output_path
            summary_rows.append(
                {
                    "Patient ID": patient_id,
                    "Number of Valid Controls": num_controls,
                    "Disconnectome Map": output_path,
                }
            )

        if not patient_avg_maps:
            raise RuntimeError("Step 5 produced no patient disconnectome maps.")

        summary_csv = os.path.join(
            self.dirs["patient_avg"], "patient_disconnectome_summary.csv"
        )
        pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)
        self.logger.info(f"Patient disconnectome summary saved to: {summary_csv}")

        return patient_avg_maps

    def step6_group_calculate(self, patient_avg_maps):
        """Clinical group analysis and Discon Score calculation.

        Step 6 is required and uses --clinical_groups_csv.

        Required clinical CSV
        ---------------------
        At minimum, the CSV must contain a patient-ID column and a group column.
        Defaults are:

            patient_id,group
            patient_001,unimpaired
            patient_002,delay
            ...

        Column names and group labels can be changed via CLI arguments.

        Calculation
        -----------
        Let D_i(x) denote the Step-5 disconnectome probability map for patient i.

        1) Unimpaired group-average map:

               U(x) = mean_i[D_i(x)], i in unimpaired group

        2) Delay group-average map:

               G(x) = mean_i[D_i(x)], i in delay group

        3) Target threshold:

               T = max_x U(x)

        4) Delay-associated target:

               Target(x) = 1, if G(x) >= T
                           0, otherwise

        5) Individual patient binary disconnectome:

               P_i(x) = 1, if D_i(x) > 0
                        0, otherwise

        6) Spatial overlap:

               Overlap_i(x) = P_i(x) AND Target(x)

        7) Discon Score:

               DisconScore_i = |P_i ∩ Target| / |Target|

        Therefore, Discon Score is the fraction of the target covered by the
        patient's disconnectome. It is bounded between 0 and 1.
        """

        clinical_csv = self.config.get("clinical_groups_csv")

        self.logger.info("Discon-Step 6: Clinical group analysis...")
        self.logger.info(
            "Step 6 formula: Discon Score = "
            "overlap(patient disconnectome binary, target) / target volume"
        )

        if not clinical_csv:
            raise ValueError(
                "Step 6 requires --clinical_groups_csv to calculate the final "
                "Discon Score."
            )

        if not os.path.exists(clinical_csv):
            raise FileNotFoundError(
                f"Clinical grouping CSV not found: {clinical_csv}"
            )

        patient_id_column = self.config.get("patient_id_column", "patient_id")
        group_column = self.config.get("group_column", "group")
        unimpaired_label = str(
            self.config.get("unimpaired_label", "unimpaired")
        )
        delay_label = str(self.config.get("delay_label", "delay"))

        clinical_df = pd.read_csv(clinical_csv)

        missing_columns = [
            col
            for col in (patient_id_column, group_column)
            if col not in clinical_df.columns
        ]
        if missing_columns:
            raise ValueError(
                "Clinical grouping CSV is missing required column(s): "
                + ", ".join(missing_columns)
            )

        clinical_df = clinical_df[[patient_id_column, group_column]].copy()
        clinical_df[patient_id_column] = clinical_df[patient_id_column].astype(str)
        clinical_df[group_column] = clinical_df[group_column].astype(str)

        unimpaired_pids = clinical_df.loc[
            clinical_df[group_column] == unimpaired_label, patient_id_column
        ].tolist()
        delay_pids = clinical_df.loc[
            clinical_df[group_column] == delay_label, patient_id_column
        ].tolist()

        unimpaired_pids = [p for p in unimpaired_pids if p in patient_avg_maps]
        delay_pids = [p for p in delay_pids if p in patient_avg_maps]

        if not unimpaired_pids or not delay_pids:
            raise ValueError(
                "Step 6 requires at least one valid patient disconnectome map in "
                f"both groups. Found unimpaired={len(unimpaired_pids)}, "
                f"delay={len(delay_pids)}."
            )

        def calculate_group_average(patient_ids):
            sum_data = None
            count = 0
            ref_img = None

            for patient_id in patient_ids:
                image = nib.load(patient_avg_maps[patient_id])
                data = image.get_fdata(dtype=np.float32)

                if sum_data is None:
                    sum_data = data.copy()
                    ref_img = image
                else:
                    if data.shape != sum_data.shape:
                        raise RuntimeError(
                            f"Shape mismatch while calculating group average: "
                            f"{patient_id}"
                        )
                    sum_data += data
                count += 1

            return sum_data / float(count), ref_img

        unimpaired_avg_data, ref_nii = calculate_group_average(unimpaired_pids)
        delay_avg_data, _ = calculate_group_average(delay_pids)

        # Save both group-average maps for traceability.
        unimpaired_avg_path = os.path.join(
            self.dirs["final_output"], "Unimpaired_Group_Average.nii.gz"
        )
        delay_avg_path = os.path.join(
            self.dirs["final_output"], "Delay_Group_Average.nii.gz"
        )

        avg_header = ref_nii.header.copy()
        avg_header.set_data_dtype(np.float32)
        nib.save(
            nib.Nifti1Image(
                unimpaired_avg_data.astype(np.float32),
                ref_nii.affine,
                avg_header,
            ),
            unimpaired_avg_path,
        )
        nib.save(
            nib.Nifti1Image(
                delay_avg_data.astype(np.float32),
                ref_nii.affine,
                avg_header,
            ),
            delay_avg_path,
        )

        max_unimpaired_val = float(np.max(unimpaired_avg_data))

        target_threshold = max_unimpaired_val + 0.01
        if max_unimpaired_val <= 0:
            raise ValueError(
                "The maximum value of the unimpaired group-average map is <= 0; "
                "a valid target cannot be constructed."
            )

        self.logger.info(
            "Step 6 target threshold = max(unimpaired group average) = "
            f"{max_unimpaired_val:.6f}"
        )

        target_data = (delay_avg_data >= max_unimpaired_val).astype(np.uint8)
        target_volume = int(np.count_nonzero(target_data))

        if target_volume == 0:
            raise ValueError(
                "The delay-group target contains zero voxels at the current "
                "threshold; Discon Score cannot be calculated."
            )

        target_header = ref_nii.header.copy()
        target_header.set_data_dtype(np.uint8)
        target_map_path = os.path.join(
            self.dirs["final_output"],
            "Target_DelayGroup_Binarized.nii.gz",
        )
        nib.save(
            nib.Nifti1Image(target_data, ref_nii.affine, target_header),
            target_map_path,
        )

        results = []
        for patient_id, patient_map_path in patient_avg_maps.items():
            patient_avg_data = nib.load(patient_map_path).get_fdata(
                dtype=np.float32
            )
            patient_binary = patient_avg_data > 0
            target_binary = target_data > 0

            overlap = np.logical_and(patient_binary, target_binary)

            patient_binary_volume = int(np.count_nonzero(patient_binary))
            overlap_volume = int(np.count_nonzero(overlap))
            discon_score = overlap_volume / float(target_volume)

            results.append(
                {
                    "Patient ID": patient_id,
                    "Patient Binary Vol (vox)": patient_binary_volume,
                    "Target Vol (vox)": target_volume,
                    "Overlap with Target (vox)": overlap_volume,
                    "Discon Score": discon_score,
                }
            )

        df_results = pd.DataFrame(results)
        score_csv = os.path.join(
            self.dirs["final_output"], "discon_scores.csv"
        )
        df_results.to_csv(score_csv, index=False)
        self.logger.info(f"Step 6 Discon Scores saved to: {score_csv}")

        return df_results
