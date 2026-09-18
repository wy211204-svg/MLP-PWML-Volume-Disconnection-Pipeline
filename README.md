# Structural-Disconnection-Score-Pipeline

A Python-based command line pipeline for Structural Disconnection Mapping analysis.

This pipeline integrates widely used neuroimaging tools including FSL, ANTs, and SimpleITK to automate:

* Image registration
* Structural disconnection mapping
* Structural disconnection score calculation

---

# 1. System Requirements

The pipeline requires the following environment:

### Operating System

* Linux or macOS (recommended)
* Windows is not recommended due to FSL compatibility

### External Neuroimaging Tools

Install the following tools and ensure they are available in your `PATH`.

#### FSL (FMRIB Software Library)

Required commands:

* `fslstats`
* `probtrackx2_gpu`

Official installation guide:

https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation

#### Python

* Python >= 3.8

---

# 2. Installation

Clone this repository and install dependencies.

```bash
# Clone repository
git clone https://github.com/wy211204-svg/Structural-Disconnection-Score-Pipeline.git
cd Structural-Disconnection-Score-Pipeline

# Install python dependencies
pip install -r requirements.txt
```

---

# 3. Data Structure

The patient dataset must follow the structure below:

```text
patients_dir/
├── patient_001/
│   ├── T1.nii.gz
│   └── lesion_mask.nii.gz
│
├── patient_002/
│   ├── T1.nii.gz
│   └── lesion_mask.nii.gz
│
└── ...
```

## 3. File Description

| File | Description |
|---|---|
| T1.nii.gz | T1 images |
| lesion_mask.nii.gz | Lesion segmentation from the T1 images |

---

# 4. Running the Pipeline

Run the pipeline from the terminal:

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

---

# 5. Arguments

| Argument | Description |
|---|---|
| --patients_dir | Root directory of patient data |
| --work_dir | Output workspace directory |
| --controls_T1_to_FA_dir | Directory containing control group T1 images registered to FA/diffusion space |
| --bedpostx_dir | BedpostX results directory |
| --JHU_T1 | JHU T1 template |
| --mat_dir | Directory containing control-to-JHU transformation matrices |
| --clinical_groups_csv | CSV containing the real clinical grouping information used in Step 6 |

---

# 6. Clinical Group Analysis and Discon Score

The clinical grouping CSV should contain patient IDs and real clinical group labels, for example:

```text
patient_id,group
patient_001,unimpaired
patient_002,delay
```

Step 6 calculates the unimpaired and delay group-average probabilistic disconnection maps. The threshold was set to the maximum value of the unimpaired group-average map plus 0.01 and was applied to the delay group-average map to generate the disconnection target map.

For each patient:

```text
Structural Disconnection Score = overlap between patient binary disconnection map and binary target map/ binary disconnection target map
```

---

# 7. Output

The pipeline automatically creates subdirectories inside the specified `--work_dir`.

The final results will be exported as:

```text
disconnection_score.csv
```

This file contains:

* Patient_id
* Structural disconnection score

---

# 8. Citation

If you use this pipeline in your research, please cite the relevant neuroimaging tools:

* FSL
* ANTs
* SimpleITK

You may also cite this repository if it is publicly released.

---

# 9. License

This project is intended for research purposes only.
