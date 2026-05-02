#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <nvtx3/nvToolsExt.h>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include "LATfield2.hpp"
#include "particles/LATfield2_Particles.hpp"
#include "particles/LATfield2_perfParticles.hpp"
#include <gsl/gsl_rng.h>

#ifndef VELOCITY_DECAY
#define VELOCITY_DECAY 0
#endif

using namespace std;
using namespace LATfield2;

static void print_device_memory_attrs(int device)
{
    int pageable_access = 0;
    int pageable_uses_host_pt = 0;
    int managed_memory = 0;
    int concurrent_managed = 0;
    int unified_addressing = 0;

    cudaDeviceGetAttribute(&pageable_access, cudaDevAttrPageableMemoryAccess, device);
    cudaDeviceGetAttribute(&pageable_uses_host_pt, cudaDevAttrPageableMemoryAccessUsesHostPageTables, device);
    cudaDeviceGetAttribute(&managed_memory, cudaDevAttrManagedMemory, device);
    cudaDeviceGetAttribute(&concurrent_managed, cudaDevAttrConcurrentManagedAccess, device);
    cudaDeviceGetAttribute(&unified_addressing, cudaDevAttrUnifiedAddressing, device);

    std::cerr << "proc#" << parallel.rank() << ": Device " << device
              << " attrs: pageableMemoryAccess=" << pageable_access
              << ", pageableUsesHostPageTables=" << pageable_uses_host_pt
              << ", managedMemory=" << managed_memory
              << ", concurrentManagedAccess=" << concurrent_managed
              << ", unifiedAddressing=" << unified_addressing
              << endl;
}

// create num_pcls particles in random positions
void initialize_particles(Particles<part_simple, part_simple_info, part_simple_dataType> * old_particles, perfParticles<part_simple, part_simple_info> *new_particles, uint64_t num_pcls, Real boxSize[3])
{
    // create a random number generator using GSL
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

    // initialize with fixed seed
    gsl_rng_set(rng, 1234);

    // create a particle
    part_simple pcl;

    // loop
    for (uint64_t i = 0; i < num_pcls; i++)
    {
        // set the position of the particle
        pcl.pos[0] = gsl_rng_uniform(rng) * boxSize[0];
        if (pcl.pos[0] >= boxSize[0])
            pcl.pos[0] -= boxSize[0];
        pcl.pos[1] = gsl_rng_uniform(rng) * boxSize[1];
        if (pcl.pos[1] >= boxSize[1])
            pcl.pos[1] -= boxSize[1];
        pcl.pos[2] = gsl_rng_uniform(rng) * boxSize[2];
        if (pcl.pos[2] >= boxSize[2])
            pcl.pos[2] -= boxSize[2];

        // set the velocity of the particle between -0.5 and 0.5
        pcl.vel[0] = (gsl_rng_uniform(rng) - 0.5) * 0.05;
        pcl.vel[1] = (gsl_rng_uniform(rng) - 0.5) * 0.05;
        pcl.vel[2] = (gsl_rng_uniform(rng) - 0.5) * 0.05;

        // set the ID of the particle
        pcl.ID = i;

        // add the particle to the particle handler
        old_particles->addParticle_global(pcl);
        new_particles->addParticle_global(pcl);
    }

    // free the random number generator
    gsl_rng_free(rng);

    // compact rows in the new particle handler
    //new_particles.compactRows(8);

    new_particles->updateRowBuffers();
}

// simple kick function

__host__ __device__ Real kick_function(double dt, double dx, part_simple* pcl, double* ref_dist, part_simple_info info, Field<Real> * fields[], Site * sites, int nf, double* param, double* out, int nout)
{
    Real v2 = 0.0;

    if (fields == NULL || nf == 0 || sites == NULL || ref_dist == NULL)
        return pcl->vel[0]*pcl->vel[0] + pcl->vel[1]*pcl->vel[1] + pcl->vel[2]*pcl->vel[2];

    for (int l = 0; l < 3; l++)
    {
        pcl->vel[l] -= dt * (VELOCITY_DECAY*pcl->vel[l] + ((1.-ref_dist[l]) * ((*fields[0])(sites[0]+l) - (*fields[0])(sites[0]-l)) + ref_dist[l] * ((*fields[0])(sites[0]+l+l) - (*fields[0])(sites[0]))) / (2.0*dx));
        v2 += pcl->vel[l] * pcl->vel[l];
    }

    return v2;
}

struct kick_function_struct
{
    __host__ __device__ Real operator()(double dt, double dx, part_simple* pcl, double* ref_dist, part_simple_info info, Field<Real> * fields[], Site * sites, int nf, double* param, double* out, int nout)
    {
        return kick_function(dt, dx, pcl, ref_dist, info, fields, sites, nf, param, out, nout);
    }
};

// main function

int main(int argc, char **argv)
{
    int n = 2, m = 2;
    int Ngrid[3] = {32, 32, 32};
    uint64_t Npcl = 32768;
    int benchmark_iterations = 1;
    // empty string
    string ofilename = "";

    for (int i = 0; i < argc; i++)
    {
        if (strcmp(argv[i], "-n") == 0)
        {
            n = atoi(argv[++i]);
        }
        if (strcmp(argv[i], "-m") == 0)
        {
            m = atoi(argv[++i]);
        }
        if (strcmp(argv[i], "-Ngrid") == 0)
        {
            Ngrid[0] = (Ngrid[1] = (Ngrid[2] = atoi(argv[++i])));
        }
        if (strcmp(argv[i], "-Npcl") == 0)
        {
            Npcl = (uint64_t) atoi(argv[++i]);
        }
        if (strcmp(argv[i], "-bench") == 0)
        {
            benchmark_iterations = atoi(argv[++i]);
        }
        if (strcmp(argv[i], "-o") == 0)
        {
            ofilename = argv[++i];
        }
    }

    parallel.initialize(n, m);

    int deviceCount;
	cudaGetDeviceCount(&deviceCount);

    if (deviceCount == 0)
    {
        std::cerr << "proc#" << parallel.rank() << ": No CUDA devices found" << endl;
        return 1;
    }

    for (int device = 0; device < deviceCount; ++device)
    {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, device);
        std::cerr << "proc#" << parallel.rank() << ": Device " << device << ": " << deviceProp.name << " with " << deviceProp.multiProcessorCount << " SMs, CC " << deviceProp.major << "." << deviceProp.minor << ", global memory " << deviceProp.totalGlobalMem / (1024*1024) << " MB" << endl << endl;
        print_device_memory_attrs(device);
    }

    nvtxMarkA("parsing finished");

    // create a lattice with size Ngrid^3
    Lattice * lat;

    cudaMallocManaged(&lat, sizeof(Lattice));
    new (lat) Lattice(3, Ngrid, 2);

    // create two particle handlers using the two implementations
    Particles<part_simple, part_simple_info, part_simple_dataType> particles_old;
    perfParticles<part_simple, part_simple_info> * particles_new;
    part_simple_info pcl_info;
    part_simple_dataType pcl_dataType;
    Real boxSize[3] = {1.0, 1.0, 1.0};

    cudaMallocManaged(&particles_new, sizeof(perfParticles<part_simple, part_simple_info>));
    new (particles_new) perfParticles<part_simple, part_simple_info>();

    // create two fields for the CIC projection
    Field<Real> density_old;
    Field<Real> density_new;

    // create a field for an external force potential
    Field<Real> * potential;

    cudaMallocManaged(&potential, sizeof(Field<Real>));
    new (potential) Field<Real>();

    // initialize the fields and particle handlers
    density_old.initialize(*lat, 1);
    density_new.initialize(*lat, 1);
    potential->initialize(*lat, 1);

    density_old.alloc(-1, Field<Real>::managed);
    density_new.alloc(-1, Field<Real>::managed);
    potential->alloc(-1, Field<Real>::managed);

    strcpy(pcl_info.type_name, "part_simple");
    pcl_info.mass = 1.0;
    pcl_info.relativistic = false;

    particles_old.initialize(pcl_info, pcl_dataType, lat, boxSize);
    //particles_new.initialize(pcl_info, &lat, boxSize, (uint32_t) Ngrid[0]+8, 8);
    particles_new->initialize(pcl_info, lat, boxSize, (uint64_t) (Npcl / n / m), 1024);

    COUT << "Initializing particles" << endl;

    initialize_particles(&particles_old, particles_new, Npcl, boxSize);

    nvtxMarkA("particles initialised");

    // initialize the force potential as a sperical trough

    COUT << "Initializing force potential" << endl;

    Site x(*lat);

    for (x.first(); x.test(); x.next())
    {
        int r2 = (x.coord(0) - Ngrid[0]/2)*(x.coord(0) - Ngrid[0]/2) + (x.coord(1) - Ngrid[1]/2)*(x.coord(1) - Ngrid[1]/2) + (x.coord(2) - Ngrid[2]/2)*(x.coord(2) - Ngrid[2]/2);

        if (r2 < Ngrid[0]*Ngrid[0]/4)
        {
            (*potential)(x) = -(cos(M_PI*sqrt(r2)/Ngrid[0])-1.0)*0.05;
        }
        else
        {
            (*potential)(x) = 0.0;
        }
    }

    potential->updateHalo(true);

    nvtxMarkA("potential defined");

    // perform unit tests

    COUT << "Initializations successful, starting unit tests" << endl << endl;

    COUT << "Testing the CIC projection" << endl;

    // test and benchmark the CIC projection
    COUT << " ...using the old implementation" << endl;
    projection_init(&density_old);

    for (int i = 0; i < 5; i++) // warm-up
        scalarProjectionCIC_project(&particles_old, &density_old);

    nvtxRangePushA("test of old CIC");
    parallel.barrier();
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < benchmark_iterations; i++) 
    {
        scalarProjectionCIC_project(&particles_old, &density_old);
        parallel.barrier();
    }
    auto end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    scalarProjectionCIC_comm(&density_old);
    auto benchmark_old = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    
    COUT << " ...using the new implementation" << endl;
    thrust::fill_n(thrust::device, density_new.data(), lat->sitesLocalGross(), Real(0));
    
    for (int i = 0; i < 5; i++) // warm-up
        particles_new->meshprojection_project(&density_new);

    nvtxRangePushA("test of new CIC");
    parallel.barrier();
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < benchmark_iterations; i++)
    {
        particles_new->meshprojection_project(&density_new);
        parallel.barrier();
    }
    end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    scalarProjectionCIC_comm(&density_new);
    auto benchmark_new = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    // compare the results
    Real tolerance = 1e-5;

    for (x.first(); x.test(); x.next())
    {
        if (fabs((density_old(x) + tolerance) / (density_new(x) + tolerance) - 1.0) > tolerance)
        {
            cout << "Error: CIC projection differs at site " << x.coord(0) << " " << x.coord(1) << " " << x.coord(2) << " --- Old: " << density_old(x) << " New: " << density_new(x) << " Ratio: " << density_old(x)/density_new(x) << endl;
            return 1;
        }
    }

    COUT << "CIC projection successful" << endl;
    COUT << "Benchmark: old implementation " << benchmark_old << " us, new implementation " << benchmark_new << " us, speed-up " << (float) benchmark_old / (float) benchmark_new << endl << endl;

    COUT << "Testing halo update" << endl;

    // test the halo update
    COUT << " ...using the old implementation" << endl;

    nvtxRangePushA("test of old halo update");
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < benchmark_iterations; i++)
    {
        density_old.updateHalo(true);
    }
    end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    benchmark_old = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    COUT << " ...using the new implementation" << endl;

    nvtxRangePushA("test of new halo update");
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < benchmark_iterations; i++)
    {
        density_new.updateHalo();
    }
    end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    benchmark_new = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    // compare the results
    for (x.haloFirst(); x.haloTest(); x.haloNext())
    {
        if (fabs((density_old(x) + tolerance) / (density_new(x) + tolerance) - 1.0) > tolerance)
        {
            cout << "Error: Halo update differs at site " << x.coord(0) << " " << x.coord(1) << " " << x.coord(2) << " --- Old: " << density_old(x) << " New: " << density_new(x) << " Ratio: " << density_old(x)/density_new(x) << endl;
            return 1;
        }
    }

    COUT << "Halo update successful" << endl;
    COUT << "Benchmark: old implementation " << benchmark_old << " us, new implementation " << benchmark_new << " us, speed-up " << (float) benchmark_old / (float) benchmark_new << endl << endl;

    COUT << "Testing particle drift" << endl;

    // test the particle drift
    COUT << " ...using the old implementation" << endl;

    nvtxRangePushA("test of old drift implementation");
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < benchmark_iterations; i++)
        particles_old.moveParticles([](double dt, double dx, part_simple* pcl,double * ref_dist, part_simple_info info, Field<Real> ** fields, Site * sites, int nf, double* param, double* out, int nout) { for (int l = 0; l < 3; l++) pcl->pos[l] += dt*pcl->vel[l]; }, 0.1);
    end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    benchmark_old = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    COUT << " ...using the new implementation" << endl;

    nvtxRangePushA("test of new drift implementation");
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < benchmark_iterations; i++)
        particles_new->moveParticles([]__host__ __device__(double dt, double dx, part_simple* pcl,double * ref_dist, part_simple_info info, Field<Real> ** fields, Site * sites, int nf, double* param, double* out, int nout) { for (int l = 0; l < 3; l++) pcl->pos[l] += dt*pcl->vel[l]; }, 0.1);
    end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    benchmark_new = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    // to check validity, do a mesh projection and compare the results
    projection_init(&density_old);
    scalarProjectionCIC_project(&particles_old, &density_old);
    scalarProjectionCIC_comm(&density_old);

    projection_init(&density_new);
    particles_new->meshprojection_project(&density_new);
    scalarProjectionCIC_comm(&density_new);

    for (x.first(); x.test(); x.next())
    {
        if (fabs((density_old(x)+tolerance)/(density_new(x)+tolerance) - 1.0) > tolerance)
        {
            cout << "Error: Particle drift differs at site " << x.coord(0) << " " << x.coord(1) << " " << x.coord(2) << " --- Old: " << density_old(x) << " New: " << density_new(x) << " Ratio: " << density_old(x)/density_new(x) << endl;
            return 1;
        }
    }

    COUT << "Particle drift successful" << endl;
    COUT << "Benchmark: old implementation " << benchmark_old << " us, new implementation " << benchmark_new << " us, speed-up " << (float) benchmark_old / (float) benchmark_new << endl << endl;

    COUT << "Testing particle kick" << endl;

    Real maxvel_old = 0.0;
    Real maxvel_new = 0.0;
    Field<Real> *potentials[1] = {potential};

    // test the particle kick
    COUT << " ...using the old implementation" << endl;
    maxvel_old = particles_old.updateVel(kick_function, 0.0);
    COUT << "Max velocity: " << maxvel_old;

    nvtxRangePushA("test of old kick implementation");
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < benchmark_iterations; i++)
    {
        maxvel_old = particles_old.updateVel(kick_function, 0.005 - (i%2)*0.01, potentials, 1);
        COUT << ", " << maxvel_old;
    }
    end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    benchmark_old = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    COUT << endl;

    COUT << " ...using the new implementation" << endl;
    maxvel_new = particles_new->updateVel(kick_function_struct(), 0.0);
    COUT << "Max velocity: " << maxvel_new;

    nvtxRangePushA("test of new kick implementation");
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < benchmark_iterations; i++)
    {
        maxvel_new = particles_new->updateVel(kick_function_struct(), 0.005 - (i%2)*0.01, potentials, 1);
        COUT << ", " << maxvel_new;
    }
    end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    benchmark_new = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    COUT << endl;

    // check if maximum velocities are within tolerance
    if (fabs(maxvel_old - maxvel_new) > tolerance)
    {
        cout << "Error: Maximum velocity differs --- Old: " << maxvel_old << " New: " << maxvel_new << " Ratio: " << maxvel_old/maxvel_new << endl;
        return 1;
    }

    COUT << "Particle kick successful" << endl;
    COUT << "Benchmark: old implementation " << benchmark_old << " us, new implementation " << benchmark_new << " us, speed-up " << (float) benchmark_old / (float) benchmark_new << endl << endl;

    COUT << "Unit tests successful" << endl << endl;

    //return 0;

    COUT << "Benchmarking the full cycle" << endl;

    // benchmark the full cycle
    COUT << " ...using the old implementation" << endl;

    nvtxRangePushA("test of old full cycle");
    start = std::chrono::high_resolution_clock::now();
    COUT << "Max velocity: " << maxvel_old;
    for (int i = 0; i < benchmark_iterations; i++)
    {
        projection_init(&density_old);
        scalarProjectionCIC_project(&particles_old, &density_old);
        scalarProjectionCIC_comm(&density_old);
        maxvel_old = particles_old.updateVel(kick_function, 0.005, potentials, 1);
        COUT << ", " << maxvel_old;
        particles_old.moveParticles([](double dt, double dx, part_simple* pcl,double * ref_dist, part_simple_info info, Field<Real> ** fields, Site * sites, int nf, double* param, double* out, int nout) { for (int l = 0; l < 3; l++) pcl->pos[l] += dt*pcl->vel[l]; }, 0.1);
    }
    end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    benchmark_old = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    COUT << endl;

    COUT << " ...using the new implementation" << endl;

    nvtxRangePushA("test of new full cycle");
    start = std::chrono::high_resolution_clock::now();
    COUT << "Max velocity: " << maxvel_new;
    for (int i = 0; i < benchmark_iterations; i++)
    {
        thrust::fill_n(thrust::device, density_new.data(), lat->sitesLocalGross(), Real(0));
        particles_new->meshprojection_project(&density_new);
        scalarProjectionCIC_comm(&density_new);
        maxvel_new = particles_new->updateVel(kick_function_struct(), 0.005, potentials, 1);
        COUT << ", " << maxvel_new;
        particles_new->moveParticles([]__host__ __device__(double dt, double dx, part_simple* pcl,double * ref_dist, part_simple_info info, Field<Real> ** fields, Site * sites, int nf, double* param, double* out, int nout) { for (int l = 0; l < 3; l++) pcl->pos[l] += dt*pcl->vel[l]; }, 0.1);
    }
    end = std::chrono::high_resolution_clock::now();
    nvtxRangePop();
    benchmark_new = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    COUT << endl;

    // check if maximum velocities are within tolerance
    if (fabs(maxvel_old - maxvel_new) > tolerance)
    {
        cout << "Error: Maximum velocity differs --- Old: " << maxvel_old << " New: " << maxvel_new << " Ratio: " << maxvel_old/maxvel_new << endl;
        return 1;
    }

    COUT << "Full cycle successful" << endl;
    COUT << "Benchmark: old implementation " << benchmark_old << " us, new implementation " << benchmark_new << " us, speed-up " << (float) benchmark_old / (float) benchmark_new << endl << endl;

    COUT << "All tests successful" << endl;

    particles_new->~perfParticles();
    cudaFree(particles_new);

    // write the results to a file
    if (ofilename != "")
    {
        density_old.saveHDF5(ofilename + "_density_old.h5");
        density_new.saveHDF5(ofilename + "_density_new.h5");
        potential->saveHDF5(ofilename + "_potential.h5");
    }

    potential->~Field();
    cudaFree(potential);

    lat->~Lattice();
    cudaFree(lat);

    return 0;
}

