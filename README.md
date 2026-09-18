# Structural Disconnection Mapping Pipeline

A Python command-line pipeline for lesion-based **structural disconnectome mapping**.

This version contains **no PWML volume-analysis module**. The previous PWML workflow, JHU white-matter atlas input, PWML result merging, and PWML-related command-line options have been removed.

The pipeline uses each patient's T1 image and lesion mask together with normative control diffusion/BedpostX data to generate a patient-level disconnectome probability/consistency map. If real clinical group assignments are supplied in the future, an optional Step 6 can construct a group-derived target and calculate an individual Discon Score.

---

## 1. Pipeline overview

### Step 1 — Patient lesion registration to each control space

For each patient and each healthy control:

1. Read the patient's `T1.nii.gz` and `lesion_mask.nii.gz`.
2. Register the patient T1 to the control T1 image using ANTs SyN.
3. Apply the forward transforms to the lesion mask.
4. Use **nearest-neighbor interpolation** for the lesion mask.

Output example:

```text
discon_1_masks_to_controls/
└── patient_001_lesion_to_control01.nii.gz
```

### Step 2 — Normative tractography with probtrackx2_gpu

The transformed lesion mask is used as the seed in each control's BedpostX model.

Current tractography parameters are retained from the original repository:

```text
-l
--onewaycondition
-c 0.15
-S 2000
--steplength=0.15
-P 5000
--fibthresh=0.01
--distthresh=0.0
--sampvox=0.0
--forcedir
--opd
```

The primary output is:

```text
fdt_paths.nii.gz
```

### Step 3 — Transform raw fdt_paths to JHU_T1 space

**This step now occurs before thresholding.**

Each raw `fdt_paths.nii.gz` map is transformed from the corresponding control space to the common `JHU_T1` space.

Per the current project specification, this transformation continues to use:

```python
sitk.sitkNearestNeighbor
```

Output example:

```text
discon_3_fdt_to_JHU/
└── patient_001_lesion_to_control01_JHU_T1_fdt_paths.nii.gz
```

### Step 4 — Threshold in JHU_T1 space

For each JHU-space `fdt_paths` map, only non-zero voxels are used to calculate:

\[
T = \mu_{nonzero} + 2\sigma_{nonzero}
\]

The binary disconnectome map is then:

\[
B(x)=
\begin{cases}
1, & fdt(x) \geq T\\
0, & fdt(x) < T
\end{cases}
\]

Output example:

```text
discon_4_binary_JHU/
├── patient_001_lesion_to_control01_JHU_T1_bin.nii.gz
└── threshold_statistics.csv
```

`threshold_statistics.csv` records the non-zero mean, standard deviation, threshold, and binary voxel count for every patient-control map.

### Step 5 — Patient-level disconnectome probability map

For a patient with `N` valid control-derived binary maps:

\[
D_p(x)=\frac{1}{N}\sum_{c=1}^{N}B_{p,c}(x)
\]

Therefore the final patient map ranges from 0 to 1.

For example:

- `1.0`: the voxel is present in the thresholded disconnectome for all valid controls.
- `0.8`: present in 80% of valid controls.
- `0.2`: present in 20% of valid controls.
- `0.0`: not present in any valid control-derived disconnectome.

Output example:

```text
discon_5_patient_disconnectome/
├── patient_001_disconnectome_prob.nii.gz
├── patient_002_disconnectome_prob.nii.gz
└── patient_disconnectome_summary.csv
```

For the current dataset, **Step 5 is the primary endpoint** when clinical grouping information is unavailable.

---

## 2. Step 6 — Clinical group target and Discon Score

Step 6 is fully defined in the code, but it is **not executed unless real clinical grouping data are provided**.

The pipeline does **not** divide patients according to filename order, patient order, or an arbitrary half split.

### 2.1 Required future clinical grouping file

Example:

```csv
patient_id,group
patient_001,unimpaired
patient_002,delay
patient_003,delay
patient_004,unimpaired
```

The default group labels are:

```text
unimpaired
delay
```

Both the column names and labels can be changed using command-line arguments.

### 2.2 Group-average disconnectome maps

Let `D_i(x)` denote the Step-5 patient disconnectome probability map.

For the unimpaired group:

\[
U(x)=\frac{1}{N_U}\sum_{i\in U}D_i(x)
\]

For the delay group:

\[
G(x)=\frac{1}{N_G}\sum_{i\in G}D_i(x)
\]

### 2.3 Target construction

The target threshold is defined as the maximum value of the unimpaired group-average map:

\[
T=\max_x U(x)
\]

The delay-associated target is:

\[
Target(x)=
\begin{cases}
1, & G(x)\geq T\\
0, & G(x)<T
\end{cases}
\]

The target volume is:

\[
V_{target}=|Target|
\]

### 2.4 Individual patient Discon Score

Each patient's Step-5 disconnectome probability map is binarized as:

\[
P_i(x)=
\begin{cases}
1, & D_i(x)>0\\
0, & D_i(x)=0
\end{cases}
\]

The spatial overlap between the patient disconnectome and target is:

\[
O_i=P_i\cap Target
\]

The final **Discon Score** is defined as:

\[
\boxed{
DisconScore_i=\frac{|P_i\cap Target|}{|Target|}
}
\]

In words:

> **Discon Score is the proportion of the target region covered by the patient's disconnectome.**

Its range is:

\[
0\leq DisconScore\leq1
\]

For example, if the target contains 1000 voxels and the patient's disconnectome overlaps 650 of those voxels:

\[
DisconScore=650/1000=0.65
\]

The patient's disconnectome outside the target does **not** increase the Discon Score.

### 2.5 Step-6 pseudocode

```python
# Real clinical groups are required.
unimpaired_maps = maps belonging to real unimpaired patients
delay_maps = maps belonging to real delay patients

unimpaired_avg = mean(unimpaired_maps)
delay_avg = mean(delay_maps)

threshold = max(unimpaired_avg)
target = delay_avg >= threshold

target_volume = count_voxels(target)

for patient in all_patients:
    patient_binary = patient_disconnectome > 0
    overlap = patient_binary AND target
    overlap_volume = count_voxels(overlap)

    discon_score = overlap_volume / target_volume
```

When `--clinical_groups_csv` is not supplied, the program logs the Step-6 definition and skips the calculation.

> Methodological note: if the Discon Score is later evaluated as a predictive biomarker, the target should be constructed within the training data/fold rather than using the test subject's outcome information.

---

## 3. System requirements

### Operating system

Linux is recommended because FSL/BedpostX/probtrackx2_gpu are required.

### External neuroimaging tools

The current pipeline directly calls:

- `fslstats`
- `probtrackx2_gpu`

The BedpostX results must already exist for all normative controls.

### Python

Python >= 3.8 is recommended.

Install Python dependencies with:

```bash
pip install -r requirements.txt
```

Python dependencies:

- NumPy
- pandas
- nibabel
- ANTsPy (`antspyx`)
- SimpleITK

---

## 4. Patient data structure

```text
patients_dir/
├── patient_001/
│   ├── T1.nii.gz
│   └── lesion_mask.nii.gz
├── patient_002/
│   ├── T1.nii.gz
│   └── lesion_mask.nii.gz
└── ...
```

Both files are required for a patient to be included.

---

## 5. Normative control inputs

### 5.1 Control T1 images aligned to diffusion/FA space

`--controls_T1_to_FA_dir` should contain control NIfTI files:

```text
controls_T1_to_FA_dir/
├── control01_xxx.nii.gz
├── control02_xxx.nii.gz
└── ...
```

The code retains the original repository convention that the filename prefix before the first underscore is used to locate the corresponding BedpostX directory.

Example:

```text
control01_xxx.nii.gz
```

maps to:

```text
bedpostx_dir/control01.bedpostX/
```

### 5.2 BedpostX directories

Each BedpostX directory must contain at least the files used by `probtrackx2_gpu`, including:

```text
merged*
nodif_brain_mask*
```

### 5.3 Control-to-JHU transforms

`--mat_dir` should contain SimpleITK-readable transform files using the naming convention:

```text
<control_id>_fwd_1.mat
```

For example:

```text
control01_xxx_fwd_1.mat
```

The transform must map the control-space `fdt_paths` image into the `JHU_T1` reference space.

---

## 6. Running the pipeline without clinical groups

This is the appropriate mode for the current dataset.

```bash
python main.py \
  --patients_dir /path/to/patients_data \
  --work_dir /path/to/output_workspace \
  --controls_T1_to_FA_dir /path/to/controls_T1_to_FA_images \
  --bedpostx_dir /path/to/bedpostx_data \
  --JHU_T1 /path/to/JHU_T1.nii.gz \
  --mat_dir /path/to/transforms_MAT
```

In this mode, Steps 1-5 run normally and Step 6 is skipped.

There is **no `--JHU_atlas` argument** because the PWML/fiber-atlas analysis has been removed.

---

## 7. Running Step 6 in the future

When real clinical grouping information becomes available:

```bash
python main.py \
  --patients_dir /path/to/patients_data \
  --work_dir /path/to/output_workspace \
  --controls_T1_to_FA_dir /path/to/controls_T1_to_FA_images \
  --bedpostx_dir /path/to/bedpostx_data \
  --JHU_T1 /path/to/JHU_T1.nii.gz \
  --mat_dir /path/to/transforms_MAT \
  --clinical_groups_csv /path/to/clinical_groups.csv
```

If different CSV columns or labels are used:

```bash
python main.py \
  ... \
  --clinical_groups_csv /path/to/clinical_groups.csv \
  --patient_id_column subject_id \
  --group_column outcome_group \
  --unimpaired_label normal \
  --delay_label delayed
```

---

## 8. Output structure

```text
work_dir/
├── pipeline.log
├── final_results.csv
├── discon_1_masks_to_controls/
├── discon_2_probtrackx/
├── discon_3_fdt_to_JHU/
├── discon_4_binary_JHU/
│   └── threshold_statistics.csv
├── discon_5_patient_disconnectome/
│   ├── patient_001_disconnectome_prob.nii.gz
│   └── patient_disconnectome_summary.csv
└── discon_6_clinical_outputs/
```

Without clinical grouping information, `discon_6_clinical_outputs/` remains empty and `final_results.csv` contains patient IDs plus the Step-5 disconnectome-map paths.

When Step 6 is executed, the directory also contains:

```text
discon_6_clinical_outputs/
├── Unimpaired_Group_Average.nii.gz
├── Delay_Group_Average.nii.gz
├── Target_DelayGroup_Binarized.nii.gz
└── discon_scores.csv
```

and `final_results.csv` additionally contains:

- patient binary disconnectome volume
- target volume
- overlap with target
- Discon Score

---

## 9. Main changes from the previous repository version

1. The entire PWML workflow has been removed.
2. `pwml.py` is no longer part of the pipeline.
3. `--JHU_atlas` has been removed.
4. Step 3 and Step 4 have been exchanged:
   - Step 3: raw `fdt_paths` -> JHU_T1 space.
   - Step 4: threshold in JHU_T1 space.
5. Step 3 continues to use nearest-neighbor interpolation by project specification.
6. Step 5 outputs patient-level `*_disconnectome_prob.nii.gz` maps.
7. The previous artificial half-split of patients into outcome groups has been removed.
8. Step 6 runs only with real clinical grouping data.
9. Discon Score has been changed from a patient-volume/target-volume ratio to a spatial target-coverage score:

\[
DisconScore=\frac{|Patient\ Disconnectome\cap Target|}{|Target|}
\]

---

## 10. Citation

If this pipeline is used in research, cite the relevant neuroimaging software and methodological sources, including FSL/probtrackx, ANTs/ANTsPy, and SimpleITK as appropriate.

---

## 11. License / research use

This pipeline is intended for research use. Verify registrations, transform directions, tractography outputs, and target construction before formal statistical analysis or publication.
