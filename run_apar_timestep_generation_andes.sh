#!/usr/bin/env bash
#
# Generate apar derivative NetCDF files on one Andes node.
#
# Submit with:
#
#   sbatch -A <project_id> run_apar_timestep_generation_andes.sh
#
# Optional overrides:
#
#   START_STEP=200 END_STEP=900 MAX_JOBS=32 \
#   GRID_FILE=../cbm18_dens3_0.5BS_516nx64ny.grid.nc \
#   sbatch -A <project_id> run_apar_timestep_generation_andes.sh
#
# By default this generates:
#
#   apar_data0_1.00200.nc ... apar_data0_1.00900.nc
#   apar_data1_1.00200.nc ... apar_data1_1.00900.nc
#
# Existing final .nc files are skipped.

#SBATCH -J apar-nc
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 04:00:00
#SBATCH -o apar_timestep_generation.%j.out
#SBATCH -e apar_timestep_generation.%j.err

set -u
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${SCRIPT_DIR}"

CONDA_ENV=${CONDA_ENV:-adios2-py}
PYTHON=${PYTHON:-python3}
CALC_SCRIPT=${CALC_SCRIPT:-calc_dxdz_dy.py}
GRID_FILE=${GRID_FILE:-../cbm18_dens3_0.5BS_516nx64ny.grid.nc}
OUTPUT_DIR=${OUTPUT_DIR:-.}
LOG_DIR=${LOG_DIR:-logs/apar_timestep_generation}
START_STEP=${START_STEP:-200}
END_STEP=${END_STEP:-900}
MAX_JOBS=${MAX_JOBS:-${SLURM_CPUS_ON_NODE:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}}

# Each Python process should use one core. Otherwise NumPy/SciPy BLAS threads can
# oversubscribe the node badly when many timesteps run in parallel.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}

DATASETS=(
  "apar_data0_1:../apar_data0_1.npy"
  "apar_data1_1:../apar_data1_1.npy"
)

die()
{
  echo "ERROR: $*" >&2
  exit 1
}

activate_conda()
{
  if ! command -v conda >/dev/null 2>&1; then
    die "conda is not on PATH. Load your conda/miniforge module before submitting, or add the module load command to this script."
  fi

  local conda_base
  conda_base=$(conda info --base) || die "conda info --base failed"
  # shellcheck disable=SC1091
  source "${conda_base}/etc/profile.d/conda.sh" || die "failed to initialize conda shell support"
  conda activate "${CONDA_ENV}" || die "failed to activate conda environment '${CONDA_ENV}'"
}

validate_inputs()
{
  [[ "${START_STEP}" =~ ^[0-9]+$ ]] || die "START_STEP must be an integer"
  [[ "${END_STEP}" =~ ^[0-9]+$ ]] || die "END_STEP must be an integer"
  [[ "${MAX_JOBS}" =~ ^[0-9]+$ ]] || die "MAX_JOBS must be an integer"
  (( END_STEP >= START_STEP )) || die "END_STEP must be >= START_STEP"
  (( MAX_JOBS >= 1 )) || die "MAX_JOBS must be >= 1"

  [[ -f "${CALC_SCRIPT}" ]] || die "missing calc script: ${CALC_SCRIPT}"
  [[ -f "${GRID_FILE}" ]] || die "missing grid file: ${GRID_FILE}"

  local spec input_file
  for spec in "${DATASETS[@]}"; do
    input_file=${spec#*:}
    [[ -f "${input_file}" ]] || die "missing input npy file: ${input_file}"
  done

  mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" || die "failed to create output/log directories"
}

run_one()
{
  local output_prefix=$1
  local input_file=$2
  local timestep=$3

  local step_string output_file tmp_file log_file
  step_string=$(printf "%05d" "${timestep}")
  output_file="${OUTPUT_DIR}/${output_prefix}.${step_string}.nc"
  tmp_file="${output_file}.tmp.${SLURM_JOB_ID:-local}.${BASHPID}"
  log_file="${LOG_DIR}/${output_prefix}.${step_string}.log"

  if [[ -f "${output_file}" ]]; then
    echo "SKIP ${output_file}"
    return 0
  fi

  rm -f "${tmp_file}"
  echo "RUN  ${output_file}"

  if "${PYTHON}" "${CALC_SCRIPT}" "${GRID_FILE}" "${input_file}" 1 "${timestep}" "${tmp_file}" >"${log_file}" 2>&1; then
    mv "${tmp_file}" "${output_file}"
    echo "DONE ${output_file}"
    return 0
  fi

  local status=$?
  rm -f "${tmp_file}"
  echo "FAIL ${output_file}; see ${log_file}" >&2
  return "${status}"
}

wait_for_batch()
{
  local failures_ref=$1
  local pid status label

  for i in "${!PIDS[@]}"; do
    pid=${PIDS[$i]}
    label=${LABELS[$i]}
    if wait "${pid}"; then
      :
    else
      status=$?
      echo "ERROR: ${label} failed with status ${status}" >&2
      printf -v "${failures_ref}" '%s' "$(( ${!failures_ref} + 1 ))"
    fi
  done

  PIDS=()
  LABELS=()
}

main()
{
  activate_conda
  validate_inputs

  echo "Work directory: ${SCRIPT_DIR}"
  echo "Conda env:      ${CONDA_ENV}"
  echo "Python:         $(command -v "${PYTHON}")"
  echo "Grid file:      ${GRID_FILE}"
  echo "Output dir:     ${OUTPUT_DIR}"
  echo "Log dir:        ${LOG_DIR}"
  echo "Timesteps:      ${START_STEP}..${END_STEP}"
  echo "Max jobs:       ${MAX_JOBS}"

  PIDS=()
  LABELS=()
  local failures=0
  local spec output_prefix input_file timestep

  for spec in "${DATASETS[@]}"; do
    output_prefix=${spec%%:*}
    input_file=${spec#*:}

    for (( timestep = START_STEP; timestep <= END_STEP; ++timestep )); do
      run_one "${output_prefix}" "${input_file}" "${timestep}" &
      PIDS+=("$!")
      LABELS+=("${output_prefix}.${timestep}")

      if (( ${#PIDS[@]} >= MAX_JOBS )); then
        wait_for_batch failures
      fi
    done
  done

  if (( ${#PIDS[@]} > 0 )); then
    wait_for_batch failures
  fi

  if (( failures > 0 )); then
    die "${failures} timestep jobs failed"
  fi

  echo "All requested timestep files are present or were generated successfully."
}

main "$@"
