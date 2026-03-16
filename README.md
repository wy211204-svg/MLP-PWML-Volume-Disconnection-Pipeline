# Neuro Pipeline: Disconnection & PWML CLI

This is a Python-based command-line tool for executing structural connectivity disconnection and Periventricular White Matter Lesion (PWML) analyses. It leverages robust neuroimaging tools like FSL, ANTs, and SimpleITK to automate image registration, tractography, cluster analysis, and clinical prognosis prediction using a built-in decision tree model.

## Prerequisites (环境依赖)

Before running this pipeline, ensure your system has the following installed:
1. **Linux/macOS** environment (required for FSL).
2. **FSL (FMRIB Software Library)**: Ensure `fslmaths`, `fslsplit`, `flirt`, `fnirt`, `applywarp`, and `cluster` are accessible in your system `PATH`.
3. **Probtrackx2 GPU**: The script requires `probtrackx2_gpu` for fast fiber tracking.
4. **Python 3.8+**

## Installation (安装指南)

Clone this repository and install Python dependencies:

git clone https://github.com/wy211204-svg/pipeline.git
cd pipeline
pip install -r requirements.txt

Data Structure (数据结构要求)

Your patient data root directory must be organized as follows:

patients_dir/
├── patient_001/
│   ├── T1_brain.nii.gz
│   ├── T1_for_dti.nii.gz
│   └── lesion_mask.nii.gz
├── patient_002/
│   ├── T1_brain.nii.gz
│   ├── ...

## Installation (安装指南)
You can run the pipeline from the terminal by providing the necessary reference directories and paths:
python main.py \
  --patients_dir /path/to/your/patients_data \
  --work_dir /path/to/output_workspace \
  --controls_t1_dir /path/to/control_t1_images \
  --target_fa /path/to/target_FA.nii.gz \
  --bedpostx_dir /path/to/bedpostx_data \
  --mean_fa /path/to/mean_FA.nii.gz \
  --mat_dir /path/to/transforms_MAT \
  --thr_map /path/to/final_thr_map.nii.gz \
  --jhu_atlas /path/to/JHU_atlas.nii.gz

Optional Flags:
--skip_discon: Skip the Disconnection pipeline.
--skip_pwml: Skip the PWML pipeline.

## Configuration (Configuration)
--patients_dir	患者数据根目录
--work_dir	输出工作目录
--controls_t1_dir	对照组 T1 图像目录
--target_fa	T1 配准目标 FA 文件路径
--bedpostx_dir	BedpostX 结果数据目录
--mean_fa	Mean FA 文件路径
--mat_dir	ANTs 变换矩阵存储目录
--thr_map	最终阈值图谱文件路径
--jhu_atlas	JHU 概率图谱文件路径


## Installation (安装指南)
The tool will create several subdirectories in your defined --work_dir. The final numerical results and the Decision Tree prediction will be exported as a CSV file:
final_results.csv
