# Makefile for compiling and running the unit tests

# Compiler
CXX = nvcc

# Compiler flags
CXXFLAGS = -std=c++17 -O2 -g -ccbin mpic++ -arch=sm_90 --extended-lambda -Xcompiler -fopenmp #-Wall -Wextra

# Include directories (modify these as needed)
INCLUDES = -I. -I/user-environment/linux-sles15-neoverse_v2/gcc-13.3.0/hdf5-1.14.5-iyjsbrml3dbr3l7cp65dgeclqlyfcdnn/include -I/user-environment/linux-sles15-neoverse_v2/gcc-13.3.0/gsl-2.8-pjzdxlsptkmjuvnrxif5x7ellp7rab3c/include -I/user-environment/linux-sles15-neoverse_v2/gcc-13.3.0/fftw-3.3.10-3yw4wbosrsa2257uitrgpge6a3mfw7ck/include #-I/user-environment/linux-sles15-neoverse_v2/gcc-13.2.0/hdf5-1.14.3-nc5ej2u5ldacfl65vfkjzol2yvib23jb/include -I/user-environment/linux-sles15-neoverse_v2/gcc-13.2.0/gsl-2.7.1-cl22k2td3lntsth2iii2rj5bqinau77z/include #-I/usr/local/include -I. -I/usr/include/hdf5/openmpi

# defines (compilation flags -D)
DEFINES = -DHDF5 -DH5_HAVE_PARALLEL -DVELOCITY_DECAY=2 -DSINGLE -DFFT3D

# Libraries to link against
LIBS =  -L/user-environment/linux-sles15-neoverse_v2/gcc-13.3.0/hdf5-1.14.5-iyjsbrml3dbr3l7cp65dgeclqlyfcdnn/lib -L/user-environment/linux-sles15-neoverse_v2/gcc-13.3.0/gsl-2.8-pjzdxlsptkmjuvnrxif5x7ellp7rab3c/lib -L/users/adamek/local_arm/lib -L/user-environment/linux-sles15-neoverse_v2/gcc-13.3.0/fftw-3.3.10-3yw4wbosrsa2257uitrgpge6a3mfw7ck/lib -lhdf5 -lgsl -lgslcblas -lm -lfftw3f -lcufft #-L/user-environment/linux-sles15-neoverse_v2/gcc-13.2.0/hdf5-1.14.3-nc5ej2u5ldacfl65vfkjzol2yvib23jb/lib -L/user-environment/linux-sles15-neoverse_v2/gcc-13.2.0/gsl-2.7.1-cl22k2td3lntsth2iii2rj5bqinau77z/lib -lhdf5 -lgsl -lgslcblas -lm

# Source files
SRCS = unit_tests.cu
TEST_SRCS = tests/fft_layout_test.cu

# Object files
OBJS = $(SRCS:.cu=.o)
TEST_OBJS = $(TEST_SRCS:.cu=.o)

# Executable name
TARGET = unit_tests
FFT_LAYOUT_TARGET = fft_layout_test

# Default target
all: $(TARGET) $(FFT_LAYOUT_TARGET)

# Rule for building the executables
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

$(FFT_LAYOUT_TARGET): $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) -o $(FFT_LAYOUT_TARGET) $(TEST_OBJS) $(LIBS)

# Rule for compiling source files
%.o: %.cu
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c $< -o $@

tests/%.o: tests/%.cu
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TEST_OBJS) $(TARGET) $(FFT_LAYOUT_TARGET)

# Run the unit tests
#srun -N 1 -n 256 --reservation=eurohack24 -C gpu -A hck ./$(TARGET) -n 16 -m 16 -Ngrid 512 -Npcl 134217728 -bench 16  2097152
run: $(TARGET)
	srun -N 1 -n 4 -C gpu -A sm97 --partition=debug ./$(TARGET) -n 2 -m 2 -Ngrid 128 -Npcl 2097152 -bench 8

profile: $(TARGET)
	srun -N 1 -n 64 --mem=64G --partition=debug -C gpu -A sm97 --hint=exclusive,nomultithread --cpu-bind=socket ./mps-wrapper.sh ./nsys_wrapper.sh ./$(TARGET) -n 8 -m 8 -Ngrid 512 -Npcl 134217728 -bench 16

#profile: $(TARGET)
#	srun -N 1 -n 4 -C gpu -A sm97 --partition=debug ./mps-wrapper.sh ./nsys_wrapper.sh ./$(TARGET) -n 2 -m 2 -Ngrid 128 -Npcl 2097152 -bench 8

#profile: $(TARGET)
#	srun -N 1 -n 64 --reservation=eurohack24 -C gpu -A hck ./mps-wrapper.sh ./nsys_wrapper.sh ./$(TARGET) -n 8 -m 8 -Ngrid 512 -Npcl 134217728 -bench 16
#profile: $(TARGET)
#	srun -N 1 -n 4 --reservation=eurohack24 -C gpu -A hck nsys profile -o unit_test_${SLURM_PROCID} ./$(TARGET) -n 2 -m 2 -Ngrid 128 -Npcl 2097152 -bench 16

.PHONY: all clean run fft-layout-test

fft-layout-test: $(FFT_LAYOUT_TARGET)
	bash tests/run_fft_layout_tests.sh
