#!/bin/bash
# Example mps-wrapper.sh usage:
# > srun [srun args] mps-wrapper.sh [cmd] [cmd args]
export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log

# added by me
export MPICH_MALLOC_FALLBACK=1
ulimit -s unlimited

# Launch MPS from a single rank per node
if [ $SLURM_LOCALID -eq 0 ]; then
    CUDA_VISIBLE_DEVICES=0,1,2,3 nvidia-cuda-mps-control -d
fi
# Wait for MPS to start
sleep 5
# Run the command

export LOCAL_RANK=$SLURM_LOCALID
export GPUS=(0 1 2 3)
export RANKS_PER_SOCKET=$(($(($SLURM_STEP_NUM_TASKS/$SLURM_STEP_NUM_NODES))/4))
export SOCKET_ID=$(($LOCAL_RANK / $RANKS_PER_SOCKET))
export GPU_ID=${GPUS[$SOCKET_ID]}
export NUMA_NODE=$GPU_ID

export CUDA_VISIBLE_DEVICES=$GPU_ID

# Set GPU memory limit (optional, adjust as needed)
#export CUDA_DEVICE_MAX_CONNECTIONS=32
#export CUDA_DEVICE_MAX_MEMORY_ALLOCATION=0.9  # Use 90% of GPU memory

# remove nsys profile when using without profiling
#if [ $SLURM_LOCALID -eq 0 ]; then
numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE "$@"
#else
# "$@"
#fi

# Quit MPS control daemon before exiting
if [ $SLURM_LOCALID -eq 0 ]; then
    echo quit | nvidia-cuda-mps-control
fi
