#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <mpi.h>
#include "LATfield2.hpp"

#ifdef FFT3D
#include <cuda_runtime.h>
#endif

using namespace LATfield2;

namespace
{
#ifdef FFT3D
void configure_cuda_device()
{
  int device_count = 0;
  cudaError_t status = cudaGetDeviceCount(&device_count);
  if (status != cudaSuccess || device_count == 0)
  {
    std::cerr << "Failed to query CUDA devices: " << cudaGetErrorString(status) << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  MPI_Comm local_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &local_comm);
  int local_rank = 0;
  MPI_Comm_rank(local_comm, &local_rank);
  MPI_Comm_free(&local_comm);

  const int device = local_rank % device_count;
  status = cudaSetDevice(device);
  if (status != cudaSuccess)
  {
    std::cerr << "cudaSetDevice(" << device << ") failed: " << cudaGetErrorString(status) << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
}
#else
void configure_cuda_device() {}
#endif

Real reference_value(const Site &site, int iteration, int boxSize)
{
  const Real freq = static_cast<Real>(iteration + 1);
  const Real norm = static_cast<Real>(boxSize);
  const Real arg0 = freq * (static_cast<Real>(site.coord(0)) + Real(0.25));
  const Real arg1 = freq * (static_cast<Real>(site.coord(1)) + Real(1.0)) * Real(2.0) / norm;
  const Real arg2 = (static_cast<Real>(site.coord(2) + iteration + 1)) * Real(0.5);
  return static_cast<Real>(std::sin(arg0)) + static_cast<Real>(0.5) * static_cast<Real>(std::cos(arg1)) + static_cast<Real>(0.25) * static_cast<Real>(std::sin(arg2));
}

Real reference_vector_value(const Site &site, int iteration, int component, int boxSize)
{
  const Real base = reference_value(site, iteration, boxSize);
  const Real modulation = static_cast<Real>(component + 1) * Real(0.125);
  return base * (Real(1.0) + modulation);
}
} // namespace

int main(int argc, char **argv)
{
  int proc_x = 1;
  int proc_y = 1;
  int boxSize = 64;
  int halo = 2;
  int iterations = 3;

  for (int i = 1; i < argc; i++)
  {
    if (argv[i][0] != '-')
      continue;
    if (strcmp(argv[i], "-n") == 0 && i + 1 < argc)
    {
      proc_x = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc)
    {
      proc_y = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc)
    {
      boxSize = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "-it") == 0 && i + 1 < argc)
    {
      iterations = atoi(argv[++i]);
    }
  }

  if (proc_x <= 1 || proc_y <= 1)
  {
    std::cerr << "FFT layout test requires n > 1 and m > 1 (received " << proc_x << "x" << proc_y << ")." << std::endl;
    return EXIT_FAILURE;
  }

  if (boxSize % proc_x != 0)
  {
    std::cerr << "FFT layout test requires Ngrid divisible by n. Ngrid=" << boxSize << ", n=" << proc_x << std::endl;
    return EXIT_FAILURE;
  }

  if (boxSize % proc_y != 0 || ((boxSize / proc_y) % 2) != 0)
  {
    std::cerr << "FFT layout test requires Ngrid divisible by m with Ngrid/m even. "
              << "Ngrid=" << boxSize << ", m=" << proc_y << std::endl;
    return EXIT_FAILURE;
  }

  parallel.initialize(proc_x, proc_y);
  configure_cuda_device();

  Lattice lat(3, boxSize, halo);
  Lattice latK;
  latK.initializeRealFFT(lat, 0);

  Field<Real> realField;
  Field<Imag> fourierField;
  Field<Real> vectorField;
  Field<Imag> vectorFourierField;
  realField.initialize(lat, 1);
  fourierField.initialize(latK, 1);
  vectorField.initialize(lat, 3);
  vectorFourierField.initialize(latK, 3);

  PlanFFT<Imag> fft_plan(&realField, &fourierField);
  PlanFFT<Imag> vector_fft_plan(&vectorField, &vectorFourierField);

  Site x(lat);
  double worst_scalar_error = 0.0;
  double worst_vector_error = 0.0;

  for (int iter = 0; iter < iterations; ++iter)
  {
    for (x.first(); x.test(); x.next())
    {
      realField(x) = reference_value(x, iter, boxSize);
      for (int c = 0; c < 3; ++c)
      {
        vectorField(x, c) = reference_vector_value(x, iter, c, boxSize);
      }
    }

    fft_plan.execute(FFT_FORWARD);
    fft_plan.execute(FFT_BACKWARD);
    vector_fft_plan.execute(FFT_FORWARD);
    vector_fft_plan.execute(FFT_BACKWARD);

    double local_scalar_error = 0.0;
    double local_vector_error = 0.0;
    for (x.first(); x.test(); x.next())
    {
      const Real expected = reference_value(x, iter, boxSize);
      const Real recovered = realField(x) / static_cast<Real>(lat.sites());
      const double delta = static_cast<double>(recovered - expected);
      if (std::abs(delta) > local_scalar_error)
      {
        local_scalar_error = std::abs(delta);
      }
      for (int c = 0; c < 3; ++c)
      {
        const Real vector_expected = reference_vector_value(x, iter, c, boxSize);
        const Real vector_recovered = vectorField(x, c) / static_cast<Real>(lat.sites());
        const double vector_delta = static_cast<double>(vector_recovered - vector_expected);
        if (std::abs(vector_delta) > local_vector_error)
        {
          local_vector_error = std::abs(vector_delta);
        }
      }
    }

    parallel.max(local_scalar_error);
    parallel.max(local_vector_error);
    if (local_scalar_error > worst_scalar_error)
    {
      worst_scalar_error = local_scalar_error;
    }
    if (local_vector_error > worst_vector_error)
    {
      worst_vector_error = local_vector_error;
    }

    if (parallel.isRoot())
    {
      std::cout << "iteration " << iter << " scalar max abs error " << local_scalar_error
                << ", vector max abs error " << local_vector_error << std::endl;
    }
  }

  parallel.max(worst_scalar_error);
  parallel.max(worst_vector_error);

  if (parallel.isRoot())
  {
    std::cout << "FFT layout test finished for layout " << proc_x << "x" << proc_y
              << " : worst scalar absolute error " << worst_scalar_error
              << ", worst vector absolute error " << worst_vector_error << std::endl;
  }

  const double tolerance = 1.0e-4;
  int failure = (worst_scalar_error > tolerance || worst_vector_error > tolerance) ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &failure, 1, MPI_INT, MPI_MAX, parallel.lat_world_comm());

  parallel.finalize();

  return failure;
}
