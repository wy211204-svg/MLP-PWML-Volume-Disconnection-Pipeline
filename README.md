# MLP-PWML-Volume-Disconnection-Pipeline

A Python-based command line pipeline for **Punctate White Matter Lesions (PWML)volume** and **Structural Disconnection Mapping** analysis.

This pipeline integrates widely used neuroimaging tools including
**FSL**, **ANTs**, and **SimpleITK** to automate:

-   Image registration
-   Fiber-specific PWML volume quantification
-   Structural disconnection score analysis


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
git clone https://github.com/wy211204-svg/MLP-PWML-Volume-Disconnection-Pipeline.git
cd pipeline

# Install python dependencies
pip install -r requirements.txt
```

------------------------------------------------------------------------

# 3. Data Structure

The patient dataset must follow the structure below:

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

## 3. File Description

| File | Description |
|-----|-------------|
| T1.nii.gz | T1 images |
| Lesion_mask.nii.gz | Lesion segmentation form the T1 images |
------------------------------------------------------------------------

# 4. Running the Pipeline

Run the pipeline from the terminal:

``` bash
python main.py --patients_dir /path/to/patients_data --work_dir /path/to/output_workspace controls_T1_to_FA_dir /path/to/controls_T1_to_FA_images --bedpostx_dir /path/to/bedpostx_data --JHU_T1 /path/to/JHU_T1.nii.gz --mat_dir /path/to/transforms_MAT --JHU_atlas /path/to/JHU_atlas.nii.gz
```

------------------------------------------------------------------------

## 5. Arguments

| Argument | Description |
|---------|-------------|
| --patients_dir | Root directory of patient data |
| --work_dir | Output workspace directory |
| --controls_T1_to_FA_dir | Directory containing control group T1 images which registered to FA maps |
| --bedpostx_dir | BedpostX results directory |
| --JHU_T1 | JHU T1 template |
| --mat_dir | Directory containing ANTs transformation matrices |
| --JHU_atlas | JHU white matter atlas |

------------------------------------------------------------------------

# 6. Optional Flags

  Flag            Description
  --------------- ---------------------------------
  --skip_pwml     Skip the PWML volume analysis
  --skip_discon   Skip the structural disconnection mapping analysis

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

-   Fiber-specific PWML volume statistics
-   structural disconnection score

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
