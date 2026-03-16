# Neuro Pipeline: Disconnection & PWML CLI

A Python-based command line pipeline for **structural disconnection
analysis** and **Periventricular White Matter Lesion (PWML) analysis**.

This pipeline integrates widely used neuroimaging tools including
**FSL**, **ANTs**, and **SimpleITK** to automate:

-   Image registration
-   Structural connectivity disconnection analysis
-   Probabilistic tractography
-   Lesion cluster analysis
-   Clinical prognosis prediction using a built‑in decision tree model

------------------------------------------------------------------------

# 1. System Requirements

The pipeline requires the following environment:

### Operating System

-   Linux or macOS (recommended)
-   Windows is not recommended due to FSL compatibility

### External Neuroimaging Tools

Install the following tools and ensure they are available in your
`PATH`.

#### FSL (FMRIB Software Library)

Required commands:

-   `fslmaths`
-   `fslsplit`
-   `flirt`
-   `fnirt`
-   `applywarp`
-   `cluster`

Official installation guide:

https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation

#### Probtrackx2 GPU

GPU version is required for efficient tractography:

-   `probtrackx2_gpu`

#### Python

-   Python \>= 3.8

------------------------------------------------------------------------

# 2. Installation

Clone this repository and install dependencies.

``` bash
# Clone repository
git clone https://github.com/wy211204-svg/pipeline.git
cd pipeline

# Install python dependencies
pip install -r requirements.txt
```

------------------------------------------------------------------------

# 3. Data Structure

The patient dataset must follow the structure below:

    patients_dir/
    ├── patient_001/
    │   ├── T1_brain.nii.gz
    │   ├── T1_for_dti.nii.gz
    │   └── lesion_mask.nii.gz
    │
    ├── patient_002/
    │   ├── T1_brain.nii.gz
    │   ├── T1_for_dti.nii.gz
    │   └── lesion_mask.nii.gz
    │
    └── ...

## 3. File Description

| File | Description |
|-----|-------------|
| T1_brain.nii.gz | Skull-stripped T1 image |
| T1_for_dti.nii.gz | T1 image used for DTI registration |
| lesion_mask.nii.gz | Lesion mask aligned with the T1 image |
------------------------------------------------------------------------

# 4. Running the Pipeline

Run the pipeline from the terminal:

``` bash
python main.py --patients_dir /path/to/patients_data --work_dir /path/to/output_workspace --controls_t1_dir /path/to/control_t1_images --target_fa /path/to/target_FA.nii.gz --bedpostx_dir /path/to/bedpostx_data --mean_fa /path/to/mean_FA.nii.gz --mat_dir /path/to/transforms_MAT --thr_map /path/to/final_thr_map.nii.gz --jhu_atlas /path/to/JHU_atlas.nii.gz
```

------------------------------------------------------------------------

## 5. Arguments

| Argument | Description |
|---------|-------------|
| --patients_dir | Root directory of patient data |
| --work_dir | Output workspace directory |
| --controls_t1_dir | Directory containing control group T1 images |
| --target_fa | Target FA template for registration |
| --bedpostx_dir | BedpostX results directory |
| --mean_fa | Mean FA template |
| --mat_dir | Directory containing ANTs transformation matrices |
| --thr_map | Final threshold probability map |
| --jhu_atlas | JHU white matter atlas |

------------------------------------------------------------------------

# 6. Optional Flags

  Flag            Description
  --------------- ---------------------------------
  --skip_discon   Skip the disconnection analysis
  --skip_pwml     Skip the PWML analysis

Example:

``` bash
python main.py ... --skip_pwml
```

------------------------------------------------------------------------

# 7. Output

The pipeline automatically creates subdirectories inside the specified
`--work_dir`.

The final results will be exported as:

    final_results.csv

This file contains:

-   Disconnection metrics
-   PWML statistics
-   Decision tree prediction results

------------------------------------------------------------------------

# 8. Citation

If you use this pipeline in your research, please cite the relevant
neuroimaging tools:

-   FSL
-   ANTs
-   SimpleITK

You may also cite this repository if it is publicly released.

------------------------------------------------------------------------

# 9. License

This project is intended for **research purposes only**.
