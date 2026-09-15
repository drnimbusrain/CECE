#!/bin/bash
#SBATCH --job-name=cece-earthaccess-tests
#SBATCH --account=<your-account>
#SBATCH --partition=<your-partition>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=cece-earthaccess-tests.%j.out
#SBATCH --error=cece-earthaccess-tests.%j.err

set -euo pipefail

# Run from the CECE repository root on an Ursa login node. This prepares the
# Python environment where outbound network access is available, then submits
# the Slurm job. Compute nodes should only activate the prepared environment.
#   bash scripts/ursa_earthaccess_tests.slurm
# Live Earthdata tests require ~/.netrc or EARTHDATA_USERNAME/EARTHDATA_TOKEN:
#   LIVE_EARTHDATA=1 bash scripts/ursa_earthaccess_tests.slurm

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  cd "$(dirname "${BASH_SOURCE[0]}")/.."

  module purge
  module use "${PWD}/modulefiles"
  module load cece_ursa.intel

  export CECE_EARTHACCESS_VENV="${CECE_EARTHACCESS_VENV:-${PWD}/.venv-ursa-earthaccess}"
  python -m venv "${CECE_EARTHACCESS_VENV}"
  source "${CECE_EARTHACCESS_VENV}/bin/activate"
  python -m pip install --upgrade pip setuptools wheel
  python -m pip install -e '.[cloud,test]'

  export LIVE_EARTHDATA="${LIVE_EARTHDATA:-0}"
  sbatch "$0"
  exit 0
fi

module purge
module use "${SLURM_SUBMIT_DIR}/modulefiles"
module load cece_ursa.intel

cd "${SLURM_SUBMIT_DIR}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

CECE_EARTHACCESS_VENV="${CECE_EARTHACCESS_VENV:-${SLURM_SUBMIT_DIR}/.venv-ursa-earthaccess}"
if [[ ! -x "${CECE_EARTHACCESS_VENV}/bin/python" ]]; then
  echo "Prepared Python environment not found: ${CECE_EARTHACCESS_VENV}" >&2
  echo "Run from an Ursa login node first: bash scripts/ursa_earthaccess_tests.slurm" >&2
  exit 2
fi
source "${CECE_EARTHACCESS_VENV}/bin/activate"

cmake -S . -B build-ursa-earthaccess \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_PYTHON_BINDINGS=ON \
  -DKokkos_ENABLE_SERIAL=ON
cmake --build build-ursa-earthaccess -j "${SLURM_CPUS_PER_TASK:-1}"

python - <<'PY'
import dask
import earthaccess
import fsspec
import h5netcdf
import s3fs
import xarray

print("earthaccess:", getattr(earthaccess, "__version__", "installed"))
print("xarray:", xarray.__version__)
print("h5netcdf:", h5netcdf.__version__)
print("fsspec:", fsspec.__version__)
print("s3fs:", s3fs.__version__)
print("dask:", dask.__version__)
PY

pytest_args=(tests/test_earthaccess_stream_bdsnp_megan3.py -v)
if [[ "${LIVE_EARTHDATA:-0}" == "1" ]]; then
  pytest_args+=(-m live_earthdata)
fi

python -m pytest "${pytest_args[@]}"
