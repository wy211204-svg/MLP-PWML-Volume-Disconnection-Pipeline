import os
import subprocess
import logging

class PipelineStepEmptyError(Exception):
    """当流程中的某个步骤未能产生任何有效输出时引发的自定义异常。"""
    pass

def setup_logger(work_dir):
    os.makedirs(work_dir, exist_ok=True)
    log_file = os.path.join(work_dir, "pipeline.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w', encoding='utf-8'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger()

def run_command_in_shell(cmd_string, logger):
    logger.info(f"--- SHELL EXEC: {cmd_string}")
    process = subprocess.Popen(
        cmd_string, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding='utf-8', errors='replace', env=os.environ, shell=True
    )
    output_lines = []
    for line in iter(process.stdout.readline, ''):
        stripped_line = line.strip()
        logger.debug(stripped_line)
        output_lines.append(stripped_line)
    retcode = process.wait()
    full_output = "\n".join(output_lines)
    return retcode == 0, f"Return Code: {retcode}", full_output

def get_fsl_volume(nifti_file, logger):
    if not os.path.exists(nifti_file): return 0
    try:
        cmd = ['fslstats', nifti_file, '-V']
        result = subprocess.check_output(cmd, text=True)
        return int(result.strip().split()[0])
    except (subprocess.CalledProcessError, IndexError, ValueError) as e:
        logger.warning(f"警告: 无法获取 {nifti_file} 的体积: {e}")
        return 0
