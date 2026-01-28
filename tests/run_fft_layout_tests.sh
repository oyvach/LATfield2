#!/usr/bin/env bash
set -euo pipefail

EXECUTABLE=${FFT_LAYOUT_EXEC:-./fft_layout_test}
MPIEXEC=${MPIEXEC:-srun -N 1 -C gpu -A sm97 -c 8 --time=5:00 --partition=debug}
MPIEXEC_FLAGS=${MPIEXEC_FLAGS:-}
NGRID=${FFT_LAYOUT_SIZE:-480}

if [[ ! -x "${EXECUTABLE}" ]]; then
  echo "Executable ${EXECUTABLE} not found. Build it first with 'make fft_layout_test'." >&2
  exit 1
fi

if [[ $# -gt 0 ]]; then
  LAYOUTS=("$@")
else
  LAYOUTS=("2x2" "2x4" "4x2")
fi

for layout in "${LAYOUTS[@]}"; do
  if [[ "${layout}" != *x* ]]; then
    echo "Invalid layout '${layout}'. Expected format NxM." >&2
    exit 2
  fi
  n=${layout%%x*}
  m=${layout##*x}
  if ! [[ "${n}" =~ ^[0-9]+$ && "${m}" =~ ^[0-9]+$ ]]; then
    echo "Invalid layout '${layout}'. N and M must be integers." >&2
    exit 2
  fi
  if (( n <= 1 || m <= 1 )); then
    echo "Invalid layout '${layout}': n and m must both be greater than 1." >&2
    exit 3
  fi
  if (( NGRID % n != 0 )); then
    echo "Invalid layout '${layout}': Ngrid (${NGRID}) must be divisible by n (${n})." >&2
    exit 4
  fi
  if (( NGRID % m != 0 )); then
    echo "Invalid layout '${layout}': Ngrid (${NGRID}) must be divisible by m (${m})." >&2
    exit 5
  fi
  if ((((NGRID / m) % 2) != 0)); then
    echo "Invalid layout '${layout}': Ngrid/m must be an even integer (currently $((${NGRID} / ${m})))." >&2
    exit 6
  fi
  np=$((n * m))
  echo "Running FFT layout test for layout ${n}x${m} using ${np} MPI ranks"
  ${MPIEXEC} ${MPIEXEC_FLAGS} -n ${np} "${EXECUTABLE}" -n ${n} -m ${m} -N ${NGRID}
done

echo "All FFT layout tests finished successfully."
