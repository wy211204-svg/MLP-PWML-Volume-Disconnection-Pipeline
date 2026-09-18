import logging
import os
import subprocess


class PipelineStepEmptyError(Exception):
    """Raised when a pipeline step fails to produce any valid output."""

    pass


def setup_logger(work_dir):
    os.makedirs(work_dir, exist_ok=True)
    log_file = os.path.join(work_dir, "pipeline.log")

    # Use a named logger to avoid duplicate handlers when imported interactively.
    logger = logging.getLogger("structural_disconnectome_pipeline")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if logger.handlers:
        logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s"
    )

    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


def run_command_in_shell(cmd_string, logger):
    logger.info(f"--- SHELL EXEC: {cmd_string}")

    process = subprocess.Popen(
        cmd_string,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=os.environ,
        shell=True,
    )

    output_lines = []
    for line in iter(process.stdout.readline, ""):
        stripped_line = line.strip()
        logger.debug(stripped_line)
        output_lines.append(stripped_line)

    retcode = process.wait()
    full_output = "\n".join(output_lines)
    return retcode == 0, f"Return Code: {retcode}", full_output


def get_fsl_volume(nifti_file, logger):
    """Return the number of non-zero/valid mask voxels reported by fslstats -V."""

    if not os.path.exists(nifti_file):
        return 0

    try:
        cmd = ["fslstats", nifti_file, "-V"]
        result = subprocess.check_output(cmd, text=True)
        return int(result.strip().split()[0])
    except (subprocess.CalledProcessError, IndexError, ValueError, FileNotFoundError) as exc:
        logger.warning(
            f"Unable to get the volume of {nifti_file} with fslstats: {exc}"
        )
        return 0
