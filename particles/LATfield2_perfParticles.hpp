#ifndef LATFIELD2_PERFPARTICLES_HPP
#define LATFIELD2_PERFPARTICLES_HPP

#include <algorithm>
#include <cstdlib>

#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>
#include <cub/cub.cuh>

#include <cuda/atomic>
#include <nvtx3/nvToolsExt.h>

#ifndef PERFPARTICLES_DEFAULT_EXTRACAPACITY
#define PERFPARTICLES_DEFAULT_EXTRACAPACITY 16
#endif

#define PERFPARTICLES_NGP 0
#define PERFPARTICLES_CIC 1
#define PERFPARTICLES_TSC 2

///////////////////////////////////////
// Header file for perfParticles class
///////////////////////////////////////
// This class is a performant implementation of a particle handler,
// using flat arrays instead of std::forward_list for improved
// performance and GPU portability. This class is intended to be
// used in conjunction with the LATfield2 library, and is designed
// to be as similar as possible to the existing particle handler.
//
// Features:
// - one array of particles per row of the lattice;
// - particle properties are stored in a struct of arrays, rather than an array of structs;
// - each array can be sorted according to the first coordinate of the particles;
// - particle-mesh projection loops over the arrays avoiding race conditions,
//   so that each row can be projected independently;

using namespace LATfield2;

struct tripleReal
{
    Real x;
    Real y;
    Real z;
};

// forward declarations for non-GH debug helpers
#ifndef GH
template <typename part, typename part_info>
class perfParticles;

struct RowCheckResult;
template <typename part, typename part_info>
__global__ void check_row_buffers(perfParticles<part, part_info> * pcl, RowCheckResult * res);
#endif

// cuda_realloc function to reallocate memory on device
template <typename T>
T * cuda_realloc(T * ptr, size_t old_size, size_t new_size)
{
    // Allocate new memory
    T * new_ptr;

    auto success = cudaMalloc(&new_ptr, new_size * sizeof(T));

    if (success != cudaSuccess)
    {
        std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Memory allocation failed in cuda_realloc");
    }

    // Copy the old memory to the new memory
    cudaMemcpy(new_ptr, ptr, old_size * sizeof(T), cudaMemcpyDeviceToDevice);
    // Free the old memory
    cudaFree(ptr);
    // Return the new memory
    return new_ptr;
}

// cuda_reallocManaged function to reallocate managed memory
template <typename T>
T * cuda_reallocManaged(T * ptr, size_t old_size, size_t new_size)
{
    // Allocate new memory
    T * new_ptr;

    auto success = cudaMallocManaged(&new_ptr, new_size * sizeof(T));

    if (success != cudaSuccess)
    {
        std::cerr << "cudaMallocManaged failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Memory allocation failed in cuda_reallocManaged");
    }

    // Copy the old memory to the new memory
    cudaMemcpy(new_ptr, ptr, old_size * sizeof(T), cudaMemcpyDefault);
    // Free the old memory
    cudaFree(ptr);
    // Return the new memory
    return new_ptr;
}

template <typename p_type, typename q_type, typename other_type>
struct row_buffer
{
    // Array of particle positions
    p_type *p;
    // Array of particle momenta
    q_type *q;
    // Array of other particle properties
    other_type *other;
    // Number of particles
    int count;
    // Capacity of the array
    int capacity;
    // Flag to check if the array is sorted
    bool sorted;
    // Resize the array
    void resize(int new_capacity);
    void resizeManaged(int new_capacity);
};

// resize function
template <typename p_type, typename q_type, typename other_type>
void row_buffer<p_type, q_type, other_type>::resize(int new_capacity)
{
    if (capacity < 0)
    {
        throw std::runtime_error("Trying to resize unmanaged row buffer");
    }
    // Reallocate the arrays using realloc
    p_type * new_p = (p_type *)realloc(p, new_capacity * 3 * sizeof(p_type));
    q_type * new_q = (q_type *)realloc(q, new_capacity * 3 * sizeof(q_type));
    other_type * new_other = (other_type *)realloc(other, new_capacity * sizeof(other_type));
    // Check if the reallocation was successful
    if (new_p == NULL || new_q == NULL || new_other == NULL)
    {
        // If not, throw an exception
        throw std::runtime_error("Memory allocation failed in row_buffer::resize");
    }
    // Update the pointers
    p = new_p;
    q = new_q;
    other = new_other;
    // Update capacity
    capacity = new_capacity;
}

// resizeManaged function
template <typename p_type, typename q_type, typename other_type>
void row_buffer<p_type, q_type, other_type>::resizeManaged(int new_capacity)
{
    resize(new_capacity);
}

template <typename part, typename part_info>
class perfParticles;

template <typename part, typename part_info>
__global__ void update_pointers(perfParticles<part, part_info> * pcl, unsigned long long int * send_begin);

template <typename part, typename part_info>
__global__ void compute_rows(perfParticles<part, part_info> * pcl, uint32_t * row, unsigned long long int starting_idx = 0);

__global__ void compute_row_offsets(const uint32_t * keys, unsigned long long int num_particles, unsigned long long int * offsets, uint32_t num_rows);

template <typename part, typename part_info, typename T>
__global__ void reorder_data(perfParticles<part, part_info> * pcl, unsigned long long int * indices_out, T * data_in, T * data_out, size_t stride, size_t ndata);

template <typename part, typename part_info, typename UpdateFunct>
__global__ void update_particles(perfParticles<part, part_info> * pcl, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, Real * v2, bool copyout = true);

template <typename part, typename part_info>
__global__ void project_particles(perfParticles<part, part_info> * pcl, Real * target, int projection_order, long start_idx, const long * jump, int stencil_k);

// Particle handler class
template <typename part, typename part_info>
class perfParticles
{
    public:
        // Constructor
        perfParticles();
        // Destructor
        ~perfParticles();
        // Initialize the particle handler
        void initialize(part_info part_global_info, Lattice * lat, Real boxSize[3], uint64_t initial_capacity = PERFPARTICLES_DEFAULT_EXTRACAPACITY, uint64_t extra_capacity = PERFPARTICLES_DEFAULT_EXTRACAPACITY);
        // Add a particle to the handler
        bool addParticle_global(part & newPart);
        // Update row buffers
        void updateRowBuffers(unsigned long long int * send_begin = nullptr);
        // Get the lattice
        Lattice & lattice() { return *lat_; }
        // Get the particle info structure
        part_info * parts_info() { return &part_global_info_; }
        // Get the number of particles
        uint64_t num_particles() { return num_particles_; }
        // Get the global cell coordinates for a particle
        __host__ __device__ void getPartCoord(Real * pos, int coord[3])
        {
            coord[0] = (int) floor(pos[0]*lat_size_[0]) % lat_size_[0];
            coord[1] = (int) floor(pos[1]*lat_size_[1]) % lat_size_[1];
            coord[2] = (int) floor(pos[2]*lat_size_[2]) % lat_size_[2];
        }
        __host__ __device__ void getPartCoord(part & p, int coord[3])
        {
            getPartCoord(p.pos, coord);
        }
        // Get the local cell coordinates for a particle
        __host__ __device__ void getPartCoordLocal(Real * pos, int coord[3])
        {
            getPartCoord(pos, coord);
            coord[1] -= coordSkip_[1];
            coord[2] -= coordSkip_[0];
        }
        __host__ __device__ void getPartCoordLocal(part & p, int coord[3])
        {
            getPartCoordLocal(p.pos, coord);
        }
        __host__ __device__ int computeRow(Real y, Real z)
        {
            int coord1 = static_cast<int>(floor(y*lat_size_[1]));
            int coord2 = static_cast<int>(floor(z*lat_size_[2]));

            if (coord1 < 0)
            {
                if ((y += boxSize_[1]) >= boxSize_[1])
                {
                    coord1 = static_cast<int>(floor((y-boxSize_[1])*lat_size_[1]));
                }
            }

            if (coord2 < 0)
            {
                if ((z += boxSize_[2]) >= boxSize_[2])
                {
                    coord2 = static_cast<int>(floor((z-boxSize_[2])*lat_size_[2]));
                }
            }

            coord1 -= coordSkip_[1];
            coord2 -= coordSkip_[0];

            if (coord1 < 0)
            {
                if (coord2 < 0)
                    return num_row_buffers_;
                else if (coord2 >= lat_size_local_[2])
                    return num_row_buffers_+6;
                else
                    return num_row_buffers_+3;
            }
            else if (coord1 >= lat_size_local_[1])
            {
                if (coord2 < 0)
                    return num_row_buffers_+2;
                else if (coord2 >= lat_size_local_[2])
                    return num_row_buffers_+8;
                else
                    return num_row_buffers_+5;
            }
            else if (coord2 < 0)
            {
                return num_row_buffers_+1;
            }
            else if (coord2 >= lat_size_local_[2])
            {
                return num_row_buffers_+7;
            }

            return coord1 + coord2 * lat_size_local_[1];
        }
        // Update velocities
        template <typename UpdateFunct>
        Real updateVel(UpdateFunct updateVel_funct, double dtau, Field<Real> ** fields = NULL, int nfields = 0, double * params = NULL, double * output = NULL, int * reduce_type = NULL, int noutput = 0);
        // Move particles
        template <typename UpdateFunct>
        void moveParticles(UpdateFunct move_funct, double dtau, Field<Real> ** fields = NULL, int nfields = 0, double * params = NULL, double * output = NULL, int * reduce_type = NULL, int noutput = 0);
        // density projection
        void meshprojection_project(Field<Real> * target, int projection_order = PERFPARTICLES_CIC);
        // generic particle-mesh projection
        template <typename UpdateFunct>
        void projectParticles(UpdateFunct project_funct, Field<Real> ** fields = NULL, int nfields = 0, double * params = NULL, double * output = NULL, int * reduce_type = NULL, int noutput = 0);
        template <typename UpdateFunct>
        void projectParticles_Async(UpdateFunct project_funct, Field<Real> ** fields = NULL, int nfields = 0, double * params = NULL, double * output = NULL, int * reduce_type = NULL, int noutput = 0);

        void sampleParticles(Real * pos, Real * vel, long * ID, size_t nsample)
        {
            cudaMemcpy(pos, p, nsample * 3 * sizeof(Real), cudaMemcpyDefault);
            cudaMemcpy(vel, q, nsample * 3 * sizeof(Real), cudaMemcpyDefault);
            cudaMemcpy(ID, other, nsample * sizeof(long), cudaMemcpyDefault);
        }

    protected:

        // particle info structure
        part_info part_global_info_;
        // lattice
        Lattice * lat_;
        // box size
        Real boxSize_[3];
        // grid points
        int lat_size_[3];
        // local grid points
        int lat_size_local_[3];
        // coordSkip
        int coordSkip_[2];
        // number of particles
        uint64_t num_particles_;

        // Array of particles
        row_buffer<Real, Real, long> * row_buffers_;
        // Number of row buffers
        uint32_t num_row_buffers_;
        // extra capacity
        uint32_t extra_capacity_;
        // global arrays
        Real * p;
        Real * q;
        long * other;
        uint64_t total_capacity_;
        // extra buffer for adding and communicating particles
        row_buffer<Real, Real, long> extra_buffer_[9];

        // Flag to check if rows are sorted
        bool rows_sorted_;

        void resizeGlobalBuffers(uint64_t new_total_capacity)
        {
            if (new_total_capacity == total_capacity_)
            {
                return;
            }

            if (new_total_capacity < num_particles_)
            {
                throw std::runtime_error("New capacity is less than number of particles");
            }

            p = cuda_realloc(p, total_capacity_ * 3, new_total_capacity * 3);
            q = cuda_realloc(q, total_capacity_ * 3, new_total_capacity * 3);
            other = cuda_realloc(other, total_capacity_, new_total_capacity);

            total_capacity_ = new_total_capacity;

            // update device mempool release threshold
            size_t threshold;
            cudaMemPool_t pool;

            cudaDeviceGetDefaultMemPool(&pool, 0);
            cudaMemPoolGetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &threshold);
            threshold = std::max(threshold, size_t((num_row_buffers_ + 1) * sizeof(unsigned long long int) + total_capacity_ * (5 * sizeof(Real) + 2 * sizeof(unsigned long long int))));
            cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &threshold);
        }

        void flushExtraBuffer(int idx = 4)
        {
            if (extra_buffer_[idx].count > 0)
            {
                if (num_particles_ + extra_buffer_[idx].count > total_capacity_)
                {
                    resizeGlobalBuffers(num_particles_ + extra_buffer_[idx].count + extra_capacity_);
                }

                cudaMemcpy(p + 3 * num_particles_, extra_buffer_[idx].p, extra_buffer_[idx].count * 3 * sizeof(Real), cudaMemcpyDefault);
                cudaMemcpy(q + 3 * num_particles_, extra_buffer_[idx].q, extra_buffer_[idx].count * 3 * sizeof(Real), cudaMemcpyDefault);
                cudaMemcpy(other + num_particles_, extra_buffer_[idx].other, extra_buffer_[idx].count * sizeof(long), cudaMemcpyDefault);

                num_particles_ += extra_buffer_[idx].count;

                extra_buffer_[idx].count = 0;
            }
        }

        // helper function to prepare communication
        void prepareComm(unsigned long long int * send_begin, uint32_t ** d_keys, void ** d_temp, cudaStream_t & stream);

        // helper function for sorting particles
        void computeSortIndices(uint32_t * d_keys, unsigned long long int * d_indices, void ** d_temp, cudaStream_t & stream);

        // helper function to reorder particles
        void reorderParticles(unsigned long long int * d_indices, void * d_temp, cudaStream_t & stream);

        // helper function for particle-mesh projection
        __device__ void project_particle(Real * target, int projection_order, long start_idx, const long * jump, int row, int stencil_k, int idx);

        // helper function for particle updates
        template <typename UpdateFunct>
        void updateParticles(UpdateFunct update_funct, double dtau, Field<Real> ** fields = NULL, int nfields = 0, double * params = NULL, double * output = NULL, int * reduce_type = NULL, int noutput = 0, Real * maxvel = nullptr, bool copyout = true, bool async = false);

        template <typename UpdateFunct>
        __host__ __device__ auto updateParticle(int row, int idx, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int noutput, bool copyout = true);

        // friend functions

        template <typename part2, typename part_info2>
        friend __global__ void update_pointers(perfParticles<part2, part_info2> * pcl, unsigned long long int * send_begin);

        template <typename part2, typename part_info2>
        friend __global__ void compute_rows(perfParticles<part2, part_info2> * pcl, uint32_t * row, unsigned long long int starting_idx);

        template <typename part2, typename part_info2, typename T>
        friend __global__ void reorder_data(perfParticles<part2, part_info2> * pcl, unsigned long long int * indices_out, T * data_in, T * data_out, size_t stride, size_t ndata);

        template <typename part2, typename part_info2, typename UpdateFunct>
        friend __global__ void update_particles(perfParticles<part2, part_info2> * pcl, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, Real * v2, bool copyout);

        template <typename part2, typename part_info2>
        friend __global__ void project_particles(perfParticles<part2, part_info2> * pcl, Real * target, int projection_order, long start_idx, const long * jump, int stencil_k);

    #ifndef GH
        template <typename part2, typename part_info2>
        friend __global__ void check_row_buffers(perfParticles<part2, part_info2> * pcl, RowCheckResult * res);
    #endif
};


// Constructor
template <typename part, typename part_info>
perfParticles<part, part_info>::perfParticles()
{
    // Initialize the number of particles
    num_particles_ = 0;
    // Initialize the row buffers
    row_buffers_ = nullptr;
    // Initialize the global arrays
    p = nullptr;
    q = nullptr;
    other = nullptr;
    total_capacity_ = 0;
    // Initialize the extra buffers
    for (int i = 0; i < 9; i++)
    {
        extra_buffer_[i].p = nullptr;
        extra_buffer_[i].q = nullptr;
        extra_buffer_[i].other = nullptr;
        extra_buffer_[i].count = 0;
        extra_buffer_[i].capacity = 0;
    }
}

// Destructor
template <typename part, typename part_info>
perfParticles<part, part_info>::~perfParticles()
{
    // Free the row buffers
    if (row_buffers_ != nullptr)
    {
        cudaFree(row_buffers_);
        cudaFree(p);
        cudaFree(q);
        cudaFree(other);
    }
    // Free the extra buffers
    for (int i = 0; i < 9; i++)
    {
        free(extra_buffer_[i].p);
        free(extra_buffer_[i].q);
        free(extra_buffer_[i].other);
    }
}

// Initialize the particle handler
template <typename part, typename part_info>
void perfParticles<part, part_info>::initialize(part_info part_global_info, Lattice * lat, Real boxSize[3], uint64_t initial_capacity, uint64_t extra_capacity)
{
    // Store the particle info
    part_global_info_ = part_global_info;
    // Store the lattice
    lat_ = lat;
    // Store the box size
    for (int i = 0; i < 3; i++)
    {
        boxSize_[i] = boxSize[i];
    }
    // Store the lattice size
    for (int i = 0; i < 3; i++)
    {
        lat_size_[i] = lat_->size(i);
        lat_size_local_[i] = lat_->sizeLocal(i);
    }
    // Store the coordSkip
    for (int i = 0; i < 2; i++)
    {
        coordSkip_[i] = lat_->coordSkip()[i];
    }
    // Store the extra capacity
    extra_capacity_ = extra_capacity;
    // compute number of row buffers
    num_row_buffers_ = lat_size_local_[1] * lat_size_local_[2];

    // Allocate the row buffers
    auto success = cudaMalloc(&row_buffers_, num_row_buffers_ * sizeof(row_buffer<Real, Real, long>));
    // Check if the allocation was successful
    if (success != cudaSuccess || row_buffers_ == nullptr)
    {
        std::cerr << " proc#" << parallel.rank() << " cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
        // If not, throw an exception
        throw std::runtime_error("Memory allocation failed for row_buffers_ in perfParticles::initialize");
    }
    // Initialize the row buffers
    thrust::for_each(thrust::device, row_buffers_, row_buffers_ + num_row_buffers_, [] __device__ (row_buffer<Real, Real, long> & rb)
    {
        rb.p = nullptr;
        rb.q = nullptr;
        rb.other = nullptr;
        rb.count = 0;
        rb.capacity = -1;
        rb.sorted = true;
    });

    // allocate global arrays
    total_capacity_ = initial_capacity;

    success = cudaMalloc(&p, total_capacity_ * 3 * sizeof(Real)); // this should be on device

    if (success != cudaSuccess)
    {
        std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Memory allocation failed for global position array in perfParticles::initialize");
    }

    success = cudaMalloc(&q, total_capacity_ * 3 * sizeof(Real)); // this should be on device

    if (success != cudaSuccess)
    {
        std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Memory allocation failed for global momentum array in perfParticles::initialize");
    }

    success = cudaMalloc(&other, total_capacity_ * sizeof(long)); // this should be on device

    if (success != cudaSuccess)
    {
        std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Memory allocation failed for global ID array in perfParticles::initialize");
    }

    // Initialize the number of particles
    num_particles_ = 0;
    // Initialize the flag
    rows_sorted_ = true;

    // allocate extra buffer for adding particles using malloc
    for (int i = 0; i < 9; i++)
    {
        extra_buffer_[i].count = 0;
        if (i == 4 || i % 2)
        {
            extra_buffer_[i].capacity = extra_capacity_;
        }
        else
        {
            extra_buffer_[i].capacity = extra_capacity_ / 4;
        }
        
        extra_buffer_[i].sorted = true;

        extra_buffer_[i].p = (Real *) malloc(extra_buffer_[i].capacity * 3 * sizeof(Real));
        extra_buffer_[i].q = (Real *) malloc(extra_buffer_[i].capacity * 3 * sizeof(Real));
        extra_buffer_[i].other = (long *) malloc(extra_buffer_[i].capacity * sizeof(long));
    }

    // update device mempool release threshold
    size_t threshold;
    cudaMemPool_t pool;

    cudaDeviceGetDefaultMemPool(&pool, 0);
    cudaMemPoolGetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &threshold);
    threshold = std::max(threshold, size_t((num_row_buffers_ + 1) * sizeof(unsigned long long int) + total_capacity_ * (5 * sizeof(Real) + 2 * sizeof(unsigned long long int))));
    cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &threshold);
}

// Add a particle to the handler
template <typename part, typename part_info>
bool perfParticles<part, part_info>::addParticle_global(part & newPart)
{
    // Check if particle is within domain
    Site x(*lat_);
    int coord[3];

    for (int i = 0; i < 3; i++)
    {
        if (newPart.pos[i] < 0 || newPart.pos[i] >= boxSize_[i])
        {
            std::cerr << "Particle out of domain" << std::endl;
            return false;
        }
    }

    getPartCoord(newPart, coord);

    if (!x.setCoord(coord))
    {
        // If not, return false
        return false;
    }

    // check that particle is in a valid row
    if (computeRow(newPart.pos[1], newPart.pos[2]) >= num_row_buffers_)
    {
        std::cerr << "Invalid row coordinates: " << newPart.pos[1] << ", " << newPart.pos[2] << endl;
        return false;
    }

    if (extra_buffer_[4].count >= extra_buffer_[4].capacity)
    {
        flushExtraBuffer(4);
    }

    // Add the particle to the buffer
    for (int i = 0; i < 3; i++)
    {
        extra_buffer_[4].p[3 * extra_buffer_[4].count + i] = newPart.pos[i];
        extra_buffer_[4].q[3 * extra_buffer_[4].count + i] = newPart.vel[i];
    }
    extra_buffer_[4].other[extra_buffer_[4].count] = newPart.ID;

    extra_buffer_[4].count++;

    rows_sorted_ = false;
    
    return true;
}

// binary search kernel for updating pointers into the global buffers
template <typename part, typename part_info>
__global__ void update_pointers(perfParticles<part, part_info> * pcl, unsigned long long int * send_begin)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < pcl->num_row_buffers_)
    {
        unsigned long long int left = 0;
        unsigned long long int right = pcl->num_particles_;

        while (left < right)
        {
            unsigned long long int mid = (left + right) / 2;

            if (pcl->computeRow(pcl->p[3*mid+1], pcl->p[3*mid+2]) < row)
            {
                left = mid + 1;
            }
            else
            {
                right = mid;
            }
        }

        pcl->row_buffers_[row].p = pcl->p + 3 * left;
        pcl->row_buffers_[row].q = pcl->q + 3 * left;
        pcl->row_buffers_[row].other = pcl->other + left;

        right = pcl->num_particles_;

        while (left < right)
        {
            unsigned long long int mid = (left + right) / 2;

            if (pcl->computeRow(pcl->p[3*mid+1], pcl->p[3*mid+2]) <= row)
            {
                left = mid + 1;
            }
            else
            {
                right = mid;
            }
        }

        pcl->row_buffers_[row].count = left - (pcl->row_buffers_[row].other - pcl->other);

        pcl->row_buffers_[row].sorted = true;
    }
    else if (send_begin != nullptr && row < pcl->num_row_buffers_ + 10)
    {
        unsigned long long int left = 0;
        unsigned long long int right = pcl->num_particles_;

        while (left < right)
        {
            unsigned long long int mid = (left + right) / 2;

            if (pcl->computeRow(pcl->p[3 * mid + 1], pcl->p[3 * mid + 2]) < row)
            {
                left = mid + 1;
            }
            else
            {
                right = mid;
            }
        }

        send_begin[row - pcl->num_row_buffers_] = left;
    }
}

// kernel to compute rows for radix sort
template <typename part, typename part_info>
__global__ void compute_rows(perfParticles<part, part_info> * pcl, uint32_t * row, unsigned long long int starting_idx)
{
    unsigned long long int idx = blockIdx.x * blockDim.x + threadIdx.x + starting_idx;

    if (idx < pcl->num_particles_)
    {
        row[idx] = static_cast<uint32_t>(pcl->computeRow(pcl->p[3*idx+1], pcl->p[3*idx+2]));
    }
}

// kernel to compute keys for the x-coordinate radix sort
#ifdef SINGLE
__global__ void compute_xkeys(uint32_t * keys, unsigned long long int * indices, Real * p, unsigned long long int num_particles)
#else
__global__ void compute_xkeys(uint64_t * keys, unsigned long long int * indices, Real * p, unsigned long long int num_particles)
#endif
{
    unsigned long long int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < num_particles)
    {
#ifdef SINGLE
        keys[idx] = __float_as_uint(p[3 * indices[idx]]);
#else
        keys[idx] = static_cast<uint64_t>(__double_as_longlong(p[3 * indices[idx]]));
#endif
    }
}

// kernel to compute row offsets from sorted row keys
// offsets size must be num_rows + 1, where offsets[num_rows] = num_particles
__global__ void compute_row_offsets(const uint32_t * keys, unsigned long long int num_particles, unsigned long long int * offsets, uint32_t num_rows)
{
    uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

    if (row > num_rows)
    {
        return;
    }

    if (row == num_rows)
    {
        offsets[row] = num_particles;
        return;
    }

    unsigned long long int left = 0;
    unsigned long long int right = num_particles;

    while (left < right)
    {
        unsigned long long int mid = (left + right) / 2;

        if (keys[mid] < row)
        {
            left = mid + 1;
        }
        else
        {
            right = mid;
        }
    }

    offsets[row] = left;
}

// kernel to reorder particle data after radix sort
template <typename part, typename part_info, typename T>
__global__ void reorder_data(perfParticles<part, part_info> * pcl, unsigned long long int * indices_out, T * data_in, T * data_out, size_t stride, size_t ndata)
{
    unsigned long long int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < ndata)
    {
        unsigned long long int pcl_idx = indices_out[idx / stride];

        data_out[idx] = data_in[pcl_idx * stride + (idx % stride)];
    }
}

// prepare communication
template <typename part, typename part_info>
void perfParticles<part, part_info>::prepareComm(unsigned long long int * send_begin, uint32_t ** d_keys, void ** d_temp, cudaStream_t & stream)
{
    unsigned long long int * send_begin_dev = send_begin;
#ifndef GH
    unsigned long long int * send_begin_tmp = nullptr;
    if (send_begin != nullptr)
    {
        auto err = cudaMallocAsync(&send_begin_tmp, 10 * sizeof(unsigned long long int), stream);
        if (err != cudaSuccess)
        {
            std::cerr << "CUDA malloc failed: " << cudaGetErrorString(err) << std::endl;
            throw std::runtime_error("Error allocating device buffer for send_begin in prepareComm");
        }
        send_begin_dev = send_begin_tmp;
    }
#endif
    uint32_t * d_keys_in = nullptr;
    unsigned long long int * d_indices_in = nullptr;
    unsigned long long int * d_indices_out = nullptr;
    size_t temp_storage_bytes = 0;
    int end_bit = static_cast<int>(ceil(log2(static_cast<double>(num_row_buffers_ + 10))));

    auto success = cudaMallocAsync(&d_indices_in, num_particles_ * sizeof(unsigned long long int) * 2L, stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_indices_in");
    }

    success = cudaMallocAsync(&d_keys_in, num_particles_ * sizeof(uint32_t), stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_keys_in");
    }

    d_indices_out = d_indices_in + num_particles_;

    // generate radix keys and indices
    nvtxRangePushA("prepareComm: compute radix keys and indices");

    if (*d_keys == nullptr)
    {
        success = cudaMallocAsync(d_keys, total_capacity_ * sizeof(uint32_t), stream);
        if (success != cudaSuccess)
        {
            std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
            throw std::runtime_error("Error in CUDA malloc for d_keys");
        }
    }

    compute_rows<<<(num_particles_+127)/128, 128, 0, stream>>>(this, d_keys_in);

    thrust::sequence(thrust::cuda::par.on(stream), d_indices_in, d_indices_in + num_particles_, 0);

    nvtxRangePop();

    // radix sort
    nvtxRangePushA("prepareComm: radix sort");
    // get temporary storage
    cub::DeviceRadixSort::SortPairs(nullptr, temp_storage_bytes, d_keys_in, *d_keys, d_indices_in, d_indices_out, num_particles_, 0, end_bit, stream);

    temp_storage_bytes = std::max(temp_storage_bytes, num_particles_ * 3 * sizeof(Real));

    success = cudaMallocAsync(d_temp, temp_storage_bytes, stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_temp_storage");
    }

    // sort
    cub::DeviceRadixSort::SortPairs(*d_temp, temp_storage_bytes, d_keys_in, *d_keys, d_indices_in, d_indices_out, num_particles_, 0, end_bit, stream);
    nvtxRangePop();

    cudaFreeAsync(d_keys_in, stream);

    reorderParticles(d_indices_out, *d_temp, stream);

    cudaFreeAsync(d_indices_in, stream);

    nvtxRangePushA("prepareComm: update pointers");

    update_pointers<<<(num_row_buffers_+137)/128, 128, 0, stream>>>(this, send_begin_dev);

    success = cudaStreamSynchronize(stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA kernel update_pointers");
    }

#ifndef GH
    if (send_begin_tmp != nullptr)
    {
        auto err = cudaMemcpyAsync(send_begin, send_begin_tmp, 10 * sizeof(unsigned long long int), cudaMemcpyDeviceToHost, stream);
        if (err != cudaSuccess)
        {
            std::cerr << "CUDA memcpy failed: " << cudaGetErrorString(err) << std::endl;
            throw std::runtime_error("Error copying send_begin from device in prepareComm");
        }
        cudaFreeAsync(send_begin_tmp, stream);
        cudaStreamSynchronize(stream);
    }
#endif
    nvtxRangePop();

    rows_sorted_ = false;
}

// compute sort indices
template <typename part, typename part_info>
void perfParticles<part, part_info>::computeSortIndices(uint32_t * d_keys, unsigned long long int * d_indices, void ** d_temp, cudaStream_t & stream)
{
    uint32_t * d_keys_out;
    unsigned long long int * d_indices_out;
    unsigned long long int * row_offsets;
#ifdef SINGLE
    uint32_t * d_xkeys_in;
    uint32_t * d_xkeys_out;
#else
    uint64_t * d_xkeys_in;
    uint64_t * d_xkeys_out;
#endif
    size_t temp_storage_bytes = 0;
    size_t temp_storage_bytes2 = 0;
    int end_bit = static_cast<int>(ceil(log2(static_cast<double>(num_row_buffers_))));

    auto success = cudaMallocAsync(&d_keys_out, num_particles_ * sizeof(uint32_t), stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_keys_out");
    }

    d_indices_out = d_indices + num_particles_;

    thrust::sequence(thrust::cuda::par.on(stream), d_indices, d_indices + num_particles_, 0);
    cub::DeviceRadixSort::SortPairs(nullptr, temp_storage_bytes, d_keys, d_keys_out, d_indices, d_indices_out, num_particles_, 0, end_bit, stream);
    temp_storage_bytes = std::max(temp_storage_bytes, num_particles_ * 3 * sizeof(Real));

    success = cudaMallocAsync(d_temp, temp_storage_bytes, stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_temp");
    }

    // sort
    cub::DeviceRadixSort::SortPairs(*d_temp, temp_storage_bytes, d_keys, d_keys_out, d_indices, d_indices_out, num_particles_, 0, end_bit, stream);

    success = cudaMallocAsync((void **) & row_offsets, (num_row_buffers_+1) * sizeof(unsigned long long int), stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for row_offsets");
    }

    // compute row offsets
    compute_row_offsets<<<(num_row_buffers_+1+127)/128, 128, 0, stream>>>(d_keys_out, num_particles_, row_offsets, num_row_buffers_);

#ifdef SINGLE
    d_xkeys_in = d_keys;
    d_xkeys_out = d_keys_out;
#else
    cudaFreeAsync(d_keys, stream);
    cudaFreeAsync(d_keys_out, stream);
    success = cudaMallocAsync((void **) & d_xkeys_in, num_particles_ * sizeof(uint64_t), stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_xkeys_in");
    }
    success = cudaMallocAsync((void **) & d_xkeys_out, num_particles_ * sizeof(uint64_t), stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_xkeys_out");
    }
#endif

    // generate keys for x-ordering
    compute_xkeys<<<(num_particles_+127)/128, 128, 0, stream>>>(d_xkeys_in, d_indices_out, p, num_particles_);
    // sort by x-keys
    cub::DeviceSegmentedSort::SortPairs(nullptr, temp_storage_bytes2, d_xkeys_in, d_xkeys_out, d_indices_out, d_indices, num_particles_, num_row_buffers_, row_offsets, row_offsets+1, stream);

    if (temp_storage_bytes2 > temp_storage_bytes)
    {
        cudaFreeAsync(*d_temp, stream);
        success = cudaMallocAsync(d_temp, temp_storage_bytes2, stream);
        if (success != cudaSuccess)
        {
            std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
            throw std::runtime_error("Error in CUDA malloc for d_temp");
        }
    }

    cub::DeviceSegmentedSort::SortPairs(*d_temp, temp_storage_bytes2, d_xkeys_in, d_xkeys_out, d_indices_out, d_indices, num_particles_, num_row_buffers_, row_offsets, row_offsets+1, stream);
    cudaFreeAsync(d_xkeys_in, stream);
    cudaFreeAsync(d_xkeys_out, stream);
    cudaFreeAsync(row_offsets, stream);
}

// reorder particles
template <typename part, typename part_info>
void perfParticles<part, part_info>::reorderParticles(unsigned long long int * d_indices, void * d_temp, cudaStream_t & stream)
{
    nvtxRangePushA("reorderParticles");

    // reorder ID data
    reorder_data<<<(num_particles_ + 127)/128, 128, 0, stream>>>(this, d_indices, other, (long *) d_temp, 1, num_particles_);

    // copy back to other
    auto success = cudaMemcpyAsync(other, d_temp, num_particles_ * sizeof(long), cudaMemcpyDeviceToDevice, stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA memcpy failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA memcpy for other array");
    }

    // reorder position data
    reorder_data<<<(num_particles_ * 3 + 127)/128, 128, 0, stream>>>(this, d_indices, p, (Real*) d_temp, 3, num_particles_ * 3);

    // reorder momentum data using old position array as target
    reorder_data<<<(num_particles_ * 3 + 127)/128, 128, 0, stream>>>(this, d_indices, q, p, 3, num_particles_ * 3);

    success = cudaStreamSynchronize(stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA kernel reorder_data");
    }

    // swap pointers around
    Real * q_old = q;
    q = p;
    p = q_old;

    // copy position data back from temporary array
    success = cudaMemcpyAsync(p, d_temp, num_particles_ * 3 * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA memcpy failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA memcpy for position array");
    }

    nvtxRangePop();
}

// Debug helper to validate row buffers on device (non-GH only, opt-in via ENABLE_ROW_CHECKS)
#ifndef GH
#ifdef ENABLE_ROW_CHECKS
struct RowCheckResult
{
    int row;
    int code;
    int start;
    int count;
};

template <typename part, typename part_info>
__global__ void check_row_buffers(perfParticles<part, part_info> * pcl, RowCheckResult * res)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;

    if (pcl == nullptr || res == nullptr)
    {
        return;
    }

    if (res->code != 0)
    {
        return;
    }

    if (pcl->row_buffers_ == nullptr)
    {
        res->row = -1;
        res->code = 4; // row_buffers_ null
        res->start = 0;
        res->count = 0;
        return;
    }

    if (row >= pcl->num_row_buffers_)
    {
        return;
    }

    auto rb = pcl->row_buffers_[row];

    if (rb.p == nullptr || rb.q == nullptr || rb.other == nullptr)
    {
        res->row = row;
        res->code = 1; // null pointer
        return;
    }

    int start = static_cast<int>((rb.p - pcl->p) / 3);
    if (start < 0 || start > pcl->num_particles_)
    {
        res->row = row;
        res->code = 2; // start out of range
        res->start = start;
        res->count = rb.count;
        return;
    }

    if (rb.count < 0 || start + rb.count > pcl->num_particles_)
    {
        res->row = row;
        res->code = 3; // count out of range
        res->start = start;
        res->count = rb.count;
        return;
    }
}
#endif
#endif

// update row buffers
template <typename part, typename part_info>
void perfParticles<part, part_info>::updateRowBuffers(unsigned long long int * send_begin)
{
    // If no particles, still refresh row buffers (non-GH) without launching device kernels
#ifndef GH
    if (num_particles_ == 0)
    {
        nvtxRangePushA("updateRowBuffers: empty fastpath");

        // populate row buffers on host and copy down
        auto * host_rb = (row_buffer<Real, Real, long> *) malloc(num_row_buffers_ * sizeof(row_buffer<Real, Real, long>));
        if (host_rb == nullptr)
        {
            throw std::runtime_error("Failed to allocate host row buffer array (empty)");
        }

        for (uint32_t r = 0; r < num_row_buffers_; ++r)
        {
            host_rb[r].p = p;
            host_rb[r].q = q;
            host_rb[r].other = other;
            host_rb[r].count = 0;
            host_rb[r].capacity = -1;
            host_rb[r].sorted = true;
        }

        auto copy_err = cudaMemcpy(row_buffers_, host_rb, num_row_buffers_ * sizeof(row_buffer<Real, Real, long>), cudaMemcpyHostToDevice);
        free(host_rb);
        if (copy_err != cudaSuccess)
        {
            std::cerr << "cudaMemcpy failed: " << cudaGetErrorString(copy_err) << std::endl;
            throw std::runtime_error("Error copying row buffers for empty case");
        }

        if (send_begin != nullptr)
        {
            unsigned long long zeros[10] = {0};
            std::copy(zeros, zeros + 10, send_begin);
        }

        rows_sorted_ = true;
        nvtxRangePop();
        return;
    }
#endif

    unsigned long long int * send_begin_tmp = nullptr;
    unsigned long long int * send_begin_dev = send_begin;
#ifndef GH
    if (send_begin != nullptr)
    {
        auto err = cudaMalloc(&send_begin_tmp, 10 * sizeof(unsigned long long int));
        if (err != cudaSuccess)
        {
            std::cerr << "cudaMalloc failed: " << cudaGetErrorString(err) << std::endl;
            throw std::runtime_error("Error allocating device buffer for send_begin");
        }
        send_begin_dev = send_begin_tmp; // kernel will fill device buffer
    }
#endif

    // check if adding buffer is empty, otherwise add particles to global buffer
    nvtxRangePushA("updateRowBuffers: add particles");

    flushExtraBuffer(4);

    nvtxRangePop();

    nvtxRangePushA("updateRowBuffers: sort");
#ifdef GH
    uint32_t * d_keys_in = nullptr;
    unsigned long long int * d_indices = nullptr;
    void * d_temp = nullptr;
    cudaStream_t pcl_stream;

    auto success = cudaStreamCreateWithFlags(&pcl_stream, cudaStreamNonBlocking);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA stream creation failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA stream creation for moveParticles");
    }

    success = cudaMallocAsync(&d_keys_in, num_particles_ * sizeof(uint32_t), pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_keys_in");
    }

    success = cudaMallocAsync(&d_indices, num_particles_ * sizeof(unsigned long long int) * 2L, pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_indices");
    }

    compute_rows<<<(num_particles_+127)/128, 128, 0, pcl_stream>>>(this, d_keys_in);

    computeSortIndices(d_keys_in, d_indices, &d_temp, pcl_stream);

    reorderParticles(d_indices, d_temp, pcl_stream);

    cudaFreeAsync(d_temp, pcl_stream);
    cudaFreeAsync(d_indices, pcl_stream);

    success = cudaStreamSynchronize(pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA kernel reorder_data");
    }

    cudaStreamDestroy(pcl_stream);
#else
#ifdef DEBUG_RADIX_SORT
    // Debug path only for non-GH: deterministic sort-by-key on rows, then reorder
    uint32_t * d_keys_in = nullptr;
    unsigned long long int * d_indices = nullptr;
    void * d_temp = nullptr;
    cudaStream_t dbg_stream;

    auto success = cudaStreamCreateWithFlags(&dbg_stream, cudaStreamNonBlocking);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA stream creation failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA stream creation for DEBUG_RADIX_SORT");
    }

    success = cudaMallocAsync(&d_keys_in, num_particles_ * sizeof(uint32_t), dbg_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_keys_in (DEBUG_RADIX_SORT)");
    }

    success = cudaMallocAsync(&d_indices, num_particles_ * sizeof(unsigned long long int), dbg_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_indices (DEBUG_RADIX_SORT)");
    }

    success = cudaMallocAsync(&d_temp, num_particles_ * 3 * sizeof(Real), dbg_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_temp (DEBUG_RADIX_SORT)");
    }

    compute_rows<<<(num_particles_+127)/128, 128, 0, dbg_stream>>>(this, d_keys_in);

    auto err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        std::cerr << "compute_rows failed: " << cudaGetErrorString(err) << std::endl;
        throw std::runtime_error("Error in compute_rows (DEBUG_RADIX_SORT)");
    }

    thrust::sequence(thrust::cuda::par.on(dbg_stream), d_indices, d_indices + num_particles_, 0);
    thrust::sort_by_key(thrust::cuda::par.on(dbg_stream), d_keys_in, d_keys_in + num_particles_, d_indices);

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        std::cerr << "thrust::sort_by_key failed: " << cudaGetErrorString(err) << std::endl;
        throw std::runtime_error("Error in thrust sort (DEBUG_RADIX_SORT)");
    }

    reorderParticles(d_indices, d_temp, dbg_stream);

    cudaFreeAsync(d_temp, dbg_stream);
    cudaFreeAsync(d_indices, dbg_stream);
    cudaFreeAsync(d_keys_in, dbg_stream);

    auto success_dbg = cudaStreamSynchronize(dbg_stream);
    if (success_dbg != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success_dbg) << std::endl;
        throw std::runtime_error("Error in DEBUG_RADIX_SORT sort/reorder path");
    }

    cudaStreamDestroy(dbg_stream);
#else
    uint32_t * d_keys_in = nullptr;
    unsigned long long int * d_indices = nullptr;
    void * d_temp = nullptr;
    cudaStream_t pcl_stream;

    auto success = cudaStreamCreateWithFlags(&pcl_stream, cudaStreamNonBlocking);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA stream creation failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA stream creation for moveParticles");
    }

    success = cudaMallocAsync(&d_keys_in, num_particles_ * sizeof(uint32_t), pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_keys_in");
    }

    success = cudaMallocAsync(&d_indices, num_particles_ * sizeof(unsigned long long int) * 2L, pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_indices");
    }

    compute_rows<<<(num_particles_+127)/128, 128, 0, pcl_stream>>>(this, d_keys_in);

    computeSortIndices(d_keys_in, d_indices, &d_temp, pcl_stream);

    reorderParticles(d_indices, d_temp, pcl_stream);

    cudaFreeAsync(d_temp, pcl_stream);
    cudaFreeAsync(d_indices, pcl_stream);

    success = cudaStreamSynchronize(pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA kernel reorder_data");
    }

    cudaStreamDestroy(pcl_stream);
#endif // DEBUG_RADIX_SORT
#endif // GH
    nvtxRangePop();

    nvtxRangePushA("updateRowBuffers: update pointers");

    update_pointers<<<(num_row_buffers_+137)/128, 128>>>(this, send_begin_dev);

    auto success2 = cudaDeviceSynchronize();

    if (success2 != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success2) << std::endl;
        throw std::runtime_error("Error in CUDA kernel update_pointers");
    }

#ifndef GH
    if (send_begin_tmp != nullptr)
    {
        auto err = cudaMemcpy(send_begin, send_begin_tmp, 10 * sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess)
        {
            std::cerr << "cudaMemcpy failed: " << cudaGetErrorString(err) << std::endl;
            throw std::runtime_error("Error copying send_begin from device");
        }
        cudaFree(send_begin_tmp);
    }
#endif

    nvtxRangePop();

    rows_sorted_ = true;
}

// update particles
template <typename part, typename part_info, typename UpdateFunct>
__global__ void update_particles(perfParticles<part, part_info> * pcl, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, Real * v2, bool copyout)
{
    int row = blockIdx.x;
    int thread_id = threadIdx.x;

    Real v2_thread = Real(0);
    double * output_thread = nullptr;
    double * output_pcl = nullptr;

    if (noutput > 0)
    {
        output_thread = (double *)alloca(noutput * sizeof(double));
        output_pcl = (double *)alloca(noutput * sizeof(double));
    }

    for (int i = 0; i < noutput; i++)
    {
        if (reduce_type[i] & (SUM | SUM_LOCAL))
        {
            output_thread[i] = 0.;
        }
        else if (reduce_type[i] & (MAX | MAX_LOCAL))
        {
            output_thread[i] = -1.e30;
        }
        else if (reduce_type[i] & (MIN | MIN_LOCAL))
        {
            output_thread[i] = 1.e30;
        }
    }

    using return_type = decltype(pcl->updateParticle(row, thread_id, update_funct, dtau, fields, nfields, params, output_thread, noutput, copyout));

    for (int idx = thread_id; idx < pcl->row_buffers_[row].count; idx += 128)
    {
        if constexpr (std::is_same_v<return_type, void>)
        {
            pcl->updateParticle(row, idx, update_funct, dtau, fields, nfields, params, output_pcl, noutput, copyout);
        }
        else
        {
            Real v2_pcl = pcl->updateParticle(row, idx, update_funct, dtau, fields, nfields, params, output_pcl, noutput, copyout);

            if (v2_pcl > v2_thread)
            {
                v2_thread = v2_pcl;
            }
        }

        for (int i = 0; i < noutput; i++)
        {
            if (reduce_type[i] & (SUM | SUM_LOCAL))
            {
                output_thread[i] +=  output_pcl[i];
            }
            else if (reduce_type[i] & (MAX | MAX_LOCAL) && output_pcl[i] > output_thread[i])
            {
                output_thread[i] = output_pcl[i];
            }
            else if (reduce_type[i] & (MIN | MIN_LOCAL) && output_pcl[i] < output_thread[i])
            {
                output_thread[i] = output_pcl[i];
            }
        }
    }

    for (int i = 0; i < noutput; i++)
    {
        cuda::atomic_ref<double, cuda::thread_scope_device> output_ref(output[i]);
        if (reduce_type[i] & SUM)
        {
            output_ref.fetch_add(output_thread[i]);
        }
        else if (reduce_type[i] & MAX)
        {
            output_ref.fetch_max(output_thread[i]);
        }
        else if (reduce_type[i] & MIN)
        {
            output_ref.fetch_min(output_thread[i]);
        }
    }

    if constexpr (!std::is_same_v<return_type, void>)
    {
        __shared__ Real smem[4];

        if (v2 != nullptr)
        {
            constexpr int W = 32;
            const int lane   = threadIdx.x & (W-1);
            const int warpId = threadIdx.x >> 5;
            constexpr int nWarps = 4; // blockDim.x = 128

            unsigned mask = __activemask();

            for (int ofs = 16; ofs > 0; ofs >>= 1)  // 16,8,4,2,1 = log2(32) steps
                v2_thread = fmax(v2_thread, __shfl_down_sync(mask, v2_thread, ofs));

            if (lane == 0)
            {
                smem[warpId] = v2_thread;
            }
            __syncthreads();

            if (warpId == 0)
            {
                v2_thread = (lane < nWarps) ? smem[lane] : Real(0);

                #pragma unroll
                for (int ofs = 2; ofs > 0; ofs >>= 1)  // 2,1 = log2(4) steps as blockDim.x = 128 => nWarps = 4
                    v2_thread = fmax(v2_thread, __shfl_down_sync(mask, v2_thread, ofs));

                if (lane == 0)
                {
                    cuda::atomic_ref<Real, cuda::thread_scope_device> v2_ref(*v2);
                    v2_ref.fetch_max(v2_thread);
                }
            }
        }
    }
}

// Particle update helper function
template <typename part, typename part_info>
template <typename UpdateFunct>
void perfParticles<part, part_info>::updateParticles(UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, Real * maxvel, bool copyout, bool async)
{
    for (int i = 0; i < noutput; i++)
    {
        if (reduce_type[i] & (SUM | SUM_LOCAL))
        {
            output[i] = 0.;
        }
        else if (reduce_type[i] & (MAX | MAX_LOCAL))
        {
            output[i] = -1.e30;
        }
        else if (reduce_type[i] & (MIN | MIN_LOCAL))
        {
            output[i] = 1.e30;
        }
    }

    Real v2 = Real(0);

    // Ensure output/reduce_type are on device; if host pointers are passed, stage through device buffers
    double * output_dev = output;
    int * reduce_type_dev = reduce_type;
    bool output_dev_owned = false;
    bool reduce_type_dev_owned = false;

    if (noutput > 0)
    {
        cudaPointerAttributes attr{};

        if (output != nullptr)
        {
            auto out_attr = cudaPointerGetAttributes(&attr, output);
            if (out_attr != cudaSuccess || attr.type != cudaMemoryTypeDevice)
            {
                auto alloc_out = cudaMalloc(&output_dev, noutput * sizeof(double));
                if (alloc_out != cudaSuccess)
                {
                    throw std::runtime_error("Failed to allocate device buffer for output");
                }
                auto copy_out = cudaMemcpy(output_dev, output, noutput * sizeof(double), cudaMemcpyHostToDevice);
                if (copy_out != cudaSuccess)
                {
                    cudaFree(output_dev);
                    throw std::runtime_error("Failed to copy output to device buffer");
                }
                output_dev_owned = true;
            }
        }

        if (reduce_type != nullptr)
        {
            auto red_attr = cudaPointerGetAttributes(&attr, reduce_type);
            if (red_attr != cudaSuccess || attr.type != cudaMemoryTypeDevice)
            {
                auto alloc_rt = cudaMalloc(&reduce_type_dev, noutput * sizeof(int));
                if (alloc_rt != cudaSuccess)
                {
                    if (output_dev_owned) cudaFree(output_dev);
                    throw std::runtime_error("Failed to allocate device buffer for reduce_type");
                }
                auto copy_rt = cudaMemcpy(reduce_type_dev, reduce_type, noutput * sizeof(int), cudaMemcpyHostToDevice);
                if (copy_rt != cudaSuccess)
                {
                    cudaFree(reduce_type_dev);
                    if (output_dev_owned) cudaFree(output_dev);
                    throw std::runtime_error("Failed to copy reduce_type to device buffer");
                }
                reduce_type_dev_owned = true;
            }
        }
    }

#ifndef GH
    // Ensure row buffers are populated before any validation or updates
    bool need_rows = !rows_sorted_;
    if (!need_rows)
    {
        if (row_buffers_ == nullptr)
        {
            need_rows = true;
        }
        else
        {
            row_buffer<Real, Real, long> rb0{};
            auto copy_rb0 = cudaMemcpy(&rb0, row_buffers_, sizeof(row_buffer<Real, Real, long>), cudaMemcpyDeviceToHost);
            if (copy_rb0 != cudaSuccess || rb0.p == nullptr || rb0.q == nullptr || rb0.other == nullptr)
            {
                need_rows = true;
            }
        }
    }

    if (need_rows)
    {
        updateRowBuffers(nullptr);

        // Re-check after updateRowBuffers
        if (row_buffers_ == nullptr)
        {
            throw std::runtime_error("row_buffers_ is null after updateRowBuffers");
        }
        row_buffer<Real, Real, long> rb0{};
        auto copy_rb0 = cudaMemcpy(&rb0, row_buffers_, sizeof(row_buffer<Real, Real, long>), cudaMemcpyDeviceToHost);
        if (copy_rb0 != cudaSuccess || rb0.p == nullptr || rb0.q == nullptr || rb0.other == nullptr)
        {
            throw std::runtime_error("Row buffer 0 is null after updateRowBuffers");
        }
    }

    // Host-validate all row buffers before launching the kernel to catch bad pointers/counts early
    {
        row_buffer<Real, Real, long> * host_rb = (row_buffer<Real, Real, long> *) malloc(num_row_buffers_ * sizeof(row_buffer<Real, Real, long>));
        if (host_rb == nullptr)
        {
            throw std::runtime_error("Failed to allocate host buffer for pre-kernel row validation");
        }

        auto copy = cudaMemcpy(host_rb, row_buffers_, num_row_buffers_ * sizeof(row_buffer<Real, Real, long>), cudaMemcpyDeviceToHost);
        if (copy != cudaSuccess)
        {
            free(host_rb);
            throw std::runtime_error("Failed to copy row buffers for pre-kernel validation");
        }

        long expected_start = 0;
        long total = 0;
        for (uint32_t r = 0; r < num_row_buffers_; ++r)
        {
            const auto & rb = host_rb[r];

            if (rb.p == nullptr || rb.q == nullptr || rb.other == nullptr)
            {
                free(host_rb);
                std::cerr << "Row buffer null pointer at row " << r << std::endl;
                throw std::runtime_error("Row buffer validation failed before update_particles (pre-kernel)");
            }

            auto start = static_cast<long>((rb.p - p) / 3);
            if (start < 0 || start > static_cast<long>(num_particles_))
            {
                free(host_rb);
                std::cerr << "Row buffer start out of range at row " << r << " start=" << start << " count=" << rb.count << std::endl;
                throw std::runtime_error("Row buffer validation failed before update_particles (pre-kernel)");
            }

            auto end = start + rb.count;
            if (rb.count < 0 || end > static_cast<long>(num_particles_))
            {
                free(host_rb);
                std::cerr << "Row buffer count out of range at row " << r << " start=" << start << " end=" << end << " num_particles=" << num_particles_ << std::endl;
                throw std::runtime_error("Row buffer validation failed before update_particles (pre-kernel)");
            }

            if (start != expected_start)
            {
                free(host_rb);
                std::cerr << "Row buffer start mismatch at row " << r << " expected_start=" << expected_start << " actual=" << start << " count=" << rb.count << std::endl;
                throw std::runtime_error("Row buffer validation failed before update_particles (pre-kernel)");
            }

            expected_start = end;
            total += rb.count;
        }

        if (total != static_cast<long>(num_particles_))
        {
            free(host_rb);
            std::cerr << "Row buffer total count mismatch total=" << total << " num_particles=" << num_particles_ << std::endl;
            throw std::runtime_error("Row buffer validation failed before update_particles (pre-kernel)");
        }

        free(host_rb);
    }
#endif

#if !defined(GH) && defined(ENABLE_ROW_CHECKS)
    // Validate row buffers before launching update_particles to catch bad pointers/counts early
    cudaPointerAttributes attr{};
    if (cudaPointerGetAttributes(&attr, row_buffers_) != cudaSuccess || attr.type != cudaMemoryTypeDevice)
    {
        throw std::runtime_error("row_buffers_ pointer is invalid or not on device");
    }
    if (cudaPointerGetAttributes(&attr, p) != cudaSuccess || attr.type != cudaMemoryTypeDevice)
    {
        throw std::runtime_error("p pointer is invalid or not on device");
    }
    if (cudaPointerGetAttributes(&attr, q) != cudaSuccess || attr.type != cudaMemoryTypeDevice)
    {
        throw std::runtime_error("q pointer is invalid or not on device");
    }
    if (cudaPointerGetAttributes(&attr, other) != cudaSuccess || attr.type != cudaMemoryTypeDevice)
    {
        throw std::runtime_error("other pointer is invalid or not on device");
    }

    // Host-side sanity check of row buffer starts/counts before launching device validator
    row_buffer<Real, Real, long> * host_rb = (row_buffer<Real, Real, long> *) malloc(num_row_buffers_ * sizeof(row_buffer<Real, Real, long>));
    if (host_rb == nullptr)
    {
        throw std::runtime_error("Failed to allocate host buffer for row validation");
    }

    auto rb_copy = cudaMemcpy(host_rb, row_buffers_, num_row_buffers_ * sizeof(row_buffer<Real, Real, long>), cudaMemcpyDeviceToHost);
    if (rb_copy != cudaSuccess)
    {
        free(host_rb);
        throw std::runtime_error("Failed to copy row_buffers_ to host for validation");
    }

    for (uint32_t r = 0; r < num_row_buffers_; ++r)
    {
        const auto & rb = host_rb[r];

        if (rb.p == nullptr || rb.q == nullptr || rb.other == nullptr)
        {
            free(host_rb);
            std::cerr << "Row buffer null pointer at row " << r << std::endl;
            throw std::runtime_error("Row buffer validation failed before update_particles (host check)");
        }

        auto start = static_cast<long>((rb.p - p) / 3);
        if (start < 0 || start > static_cast<long>(num_particles_))
        {
            free(host_rb);
            std::cerr << "Row buffer start out of range at row " << r << " start=" << start << " count=" << rb.count << std::endl;
            throw std::runtime_error("Row buffer validation failed before update_particles (host check)");
        }

        auto end = start + rb.count;
        if (rb.count < 0 || end > static_cast<long>(num_particles_))
        {
            free(host_rb);
            std::cerr << "Row buffer count out of range at row " << r << " start=" << start << " end=" << end << " num_particles=" << num_particles_ << std::endl;
            throw std::runtime_error("Row buffer validation failed before update_particles (host check)");
        }
    }

    free(host_rb);

    // Device-side validation is optional; disable by default to avoid masking real work with a debug fault
#ifdef ENABLE_ROW_CHECKS_DEVICE
    RowCheckResult * chk = nullptr;
    auto alloc_chk = cudaMalloc(&chk, sizeof(RowCheckResult));
    if (alloc_chk == cudaSuccess)
    {
        cudaMemset(chk, 0, sizeof(RowCheckResult));
        check_row_buffers<part, part_info><<<(num_row_buffers_ + 127) / 128, 128>>>(this, chk);
        auto chk_sync = cudaDeviceSynchronize();
        if (chk_sync == cudaSuccess)
        {
            RowCheckResult host_chk{};
            cudaMemcpy(&host_chk, chk, sizeof(RowCheckResult), cudaMemcpyDeviceToHost);
            if (host_chk.code != 0)
            {
                std::cerr << "Row buffer validation failed: row=" << host_chk.row
                          << " code=" << host_chk.code << " start=" << host_chk.start
                          << " count=" << host_chk.count << std::endl;
                throw std::runtime_error("Row buffer validation failed before update_particles");
            }
        }
        else
        {
            std::cerr << "check_row_buffers kernel failed: " << cudaGetErrorString(chk_sync) << std::endl;
        }
        cudaFree(chk);
    }
#endif
#endif

    update_particles<<<num_row_buffers_, 128>>>(this, update_funct, dtau, fields, nfields, params, output_dev, reduce_type_dev, noutput, &v2, copyout);

    if (async)
    {
        return;
    }

    auto success = cudaDeviceSynchronize();

    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA kernel update_particles");
    }

    // copy staged outputs back to host if needed
    if (output_dev_owned && output != nullptr)
    {
        auto copy_back = cudaMemcpy(output, output_dev, noutput * sizeof(double), cudaMemcpyDeviceToHost);
        if (copy_back != cudaSuccess)
        {
            std::cerr << "cudaMemcpy failed copying output back: " << cudaGetErrorString(copy_back) << std::endl;
            throw std::runtime_error("Failed to copy output back to host after update_particles");
        }
        cudaFree(output_dev);
    }
    if (reduce_type_dev_owned)
    {
        cudaFree(reduce_type_dev);
    }

    for (int i = 0; i < noutput; i++)
    {
        if (reduce_type[i] & SUM)
        {
            parallel.sum(output[i]);
        }
        else if (reduce_type[i] & MAX)
        {
            parallel.max(output[i]);
        }
        else if (reduce_type[i] & MIN)
        {
            parallel.min(output[i]);
        }
    }

    if (maxvel != nullptr)
    {
        *maxvel = sqrt(v2);
    }
}

// Particle update helper function
template <typename part, typename part_info>
template <typename UpdateFunct>
__host__ __device__ auto perfParticles<part, part_info>::updateParticle(int row, int idx, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int noutput, bool copyout)
{
    // Initialize sites for field operations
    Site * sites = nullptr;

    if (nfields > 0)
    {
        sites = (Site *)alloca(nfields * sizeof(Site));
    }

    Real dx = boxSize_[0] / lat_size_[0];

    int coord[3];
    double frac[3];
    part pcl;

    Real v2;

    coord[0] = (int) floor(row_buffers_[row].p[3*idx]*lat_size_[0]) % lat_size_[0];
    coord[1] = (int) floor(row_buffers_[row].p[3*idx+1]*lat_size_[1]) % lat_size_[1];
    coord[2] = (int) floor(row_buffers_[row].p[3*idx+2]*lat_size_[2]) % lat_size_[2];
        
    for (int k = 0; k < nfields; k++)
    {

        sites[k] = Site(fields[k]->lattice(), fields[k]->lattice().siteFirst() 
                        + coord[0]*fields[k]->lattice().jump(0) 
                        + (coord[1] - coordSkip_[1])*fields[k]->lattice().jump(1)
                        + (coord[2] - coordSkip_[0])*fields[k]->lattice().jump(2));
    }

    frac[0] = row_buffers_[row].p[3*idx]*lat_size_[0] - coord[0];
    frac[1] = row_buffers_[row].p[3*idx+1]*lat_size_[1] - coord[1];
    frac[2] = row_buffers_[row].p[3*idx+2]*lat_size_[2] - coord[2];

    pcl.pos[0] = row_buffers_[row].p[3*idx];
    pcl.pos[1] = row_buffers_[row].p[3*idx+1];
    pcl.pos[2] = row_buffers_[row].p[3*idx+2];
    pcl.vel[0] = row_buffers_[row].q[3*idx];
    pcl.vel[1] = row_buffers_[row].q[3*idx+1];
    pcl.vel[2] = row_buffers_[row].q[3*idx+2];
    pcl.ID = row_buffers_[row].other[idx];

    using return_type = decltype(update_funct(dtau, dx, &pcl, frac, part_global_info_, fields, sites, nfields, params, output, noutput));

    if constexpr (std::is_same_v<return_type, void>)
    {
        update_funct(dtau, dx, &pcl, frac, part_global_info_, fields, sites, nfields, params, output, noutput);
    }
    else
    {
        v2 = update_funct(dtau, dx, &pcl, frac, part_global_info_, fields, sites, nfields, params, output, noutput);
    }

    if (copyout)
    {
        row_buffers_[row].p[3*idx] = pcl.pos[0];
        row_buffers_[row].p[3*idx+1] = pcl.pos[1];
        row_buffers_[row].p[3*idx+2] = pcl.pos[2];
        row_buffers_[row].q[3*idx] = pcl.vel[0];
        row_buffers_[row].q[3*idx+1] = pcl.vel[1];
        row_buffers_[row].q[3*idx+2] = pcl.vel[2];
    }

    if constexpr (!std::is_same_v<return_type, void>)
        return v2;
}

// Project particles
template <typename part, typename part_info>
template <typename UpdateFunct>
void perfParticles<part, part_info>::projectParticles(UpdateFunct project_funct, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput)
{
    updateParticles(project_funct, 0., fields, nfields, params, output, reduce_type, noutput, nullptr, false);
}

template <typename part, typename part_info>
template <typename UpdateFunct>
void perfParticles<part, part_info>::projectParticles_Async(UpdateFunct project_funct, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput)
{
    updateParticles(project_funct, 0., fields, nfields, params, output, reduce_type, noutput, nullptr, false, true);
}

// Update velocities
template <typename part, typename part_info>
template <typename UpdateFunct>
Real perfParticles<part, part_info>::updateVel(UpdateFunct updateVel_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput)
{
    Real maxvel = 0.;

    updateParticles(updateVel_funct, dtau, fields, nfields, params, output, reduce_type, noutput, &maxvel);

    return maxvel;
}

// Move particles
template <typename part, typename part_info>
template <typename UpdateFunct>
void perfParticles<part, part_info>::moveParticles(UpdateFunct move_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput)
{
    unsigned long long int send_begin[10];
#ifndef DEBUG_RADIX_SORT
    uint32_t * d_keys = nullptr;
    void * d_temp = nullptr;
    uint64_t start_num_particles;
    uint64_t start_capacity = total_capacity_;
    cudaStream_t pcl_stream;

    auto success = cudaStreamCreateWithFlags(&pcl_stream, cudaStreamNonBlocking);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA stream creation failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA stream creation for moveParticles");
    }
#endif

    nvtxRangePushA("moveParticles: updateParticles");

    updateParticles(move_funct, dtau, fields, nfields, params, output, reduce_type, noutput);

    nvtxRangePop();

    thrust::for_each(thrust::device, row_buffers_, row_buffers_ + num_row_buffers_, [] __device__ (row_buffer<Real, Real, long> & rb) { rb.sorted = false;});

    rows_sorted_ = false;

#ifdef DEBUG_RADIX_SORT
    nvtxRangePushA("moveParticles: updateRowBuffers");

    updateRowBuffers(&send_begin[0]);

    nvtxRangePop();
#else
    nvtxRangePushA("moveParticles: prepareComm");

    prepareComm(&send_begin[0], &d_keys, &d_temp, pcl_stream);

    cudaFreeAsync(d_temp, pcl_stream); // free temporary storage for now (FIXME: maybe reuse later)

    nvtxRangePop();
#endif
    // memory movement

    nvtxRangePushA("moveParticles: set up communication buffers");

    num_particles_ = send_begin[0];
#ifndef DEBUG_RADIX_SORT
    start_num_particles = num_particles_;
#endif

    // loop over send buffers
    for (int proc = 0; proc < 9; proc++)
    {
        extra_buffer_[proc].count = send_begin[proc+1] - send_begin[proc];

        if (extra_buffer_[proc].count > 0)
        {
            // check if there is capacity in the extra buffer
            if (extra_buffer_[proc].count > extra_buffer_[proc].capacity)
            {
                if (proc == 4 || proc % 2)
                {
                    extra_buffer_[proc].resizeManaged(extra_buffer_[proc].count + extra_capacity_);
                }
                else
                {
                    extra_buffer_[proc].resizeManaged(extra_buffer_[proc].count + (extra_capacity_ / 4));
                }
            }

            // copy particles to send buffer using cudaMemcpy
            cudaMemcpy(extra_buffer_[proc].p, p + 3*send_begin[proc], extra_buffer_[proc].count * 3 * sizeof(Real), cudaMemcpyDefault);
            cudaMemcpy(extra_buffer_[proc].q, q + 3*send_begin[proc], extra_buffer_[proc].count * 3 * sizeof(Real), cudaMemcpyDefault);
            cudaMemcpy(extra_buffer_[proc].other, other + send_begin[proc], extra_buffer_[proc].count * sizeof(long), cudaMemcpyDefault);
        }
    }

    nvtxRangePop();

    // send particles to other processors

    nvtxRangePushA("moveParticles: send particles (MPI)");

    long buffer_sizes[3];

    // first, even ranks send upwards odd ranks receive
    if (parallel.grid_rank()[1] % 2 == 0) // even rank
    {
        buffer_sizes[0] = extra_buffer_[2].count;
        buffer_sizes[1] = extra_buffer_[5].count;
        buffer_sizes[2] = extra_buffer_[8].count;

        parallel.send_dim1(buffer_sizes, 3, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
        if (buffer_sizes[0] > 0)
        {
            parallel.send_dim1(extra_buffer_[2].p, 3*extra_buffer_[2].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[2].q, 3*extra_buffer_[2].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[2].other, extra_buffer_[2].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
        }
        if (buffer_sizes[1] > 0)
        {
            parallel.send_dim1(extra_buffer_[5].p, 3*extra_buffer_[5].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[5].q, 3*extra_buffer_[5].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[5].other, extra_buffer_[5].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
        }
        if (buffer_sizes[2] > 0)
        {
            parallel.send_dim1(extra_buffer_[8].p, 3*extra_buffer_[8].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[8].q, 3*extra_buffer_[8].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[8].other, extra_buffer_[8].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
        }
    }
    
    if (parallel.grid_rank()[1] % 2 == 1 || (parallel.grid_rank()[1] == 0 && parallel.grid_size()[1] % 2 == 1)) // odd rank or rank 0 with odd grid size
    {
        parallel.receive_dim1(buffer_sizes, 3, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);

        // receive particles from lower ranks
        if (buffer_sizes[0] > 0)
        {
            if (extra_buffer_[1].count + buffer_sizes[0] > extra_buffer_[1].capacity)
            {
                extra_buffer_[1].resizeManaged(extra_buffer_[1].count + buffer_sizes[0] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[1].p+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].q+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].other+extra_buffer_[1].count, buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            extra_buffer_[1].count += buffer_sizes[0];
        }
        if (buffer_sizes[1] > 0)
        {        
            if (extra_buffer_[4].count + buffer_sizes[1] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[1] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[4].p, 3*buffer_sizes[1], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[4].q, 3*buffer_sizes[1], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[4].other, buffer_sizes[1], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            extra_buffer_[4].count += buffer_sizes[1];
        }
        if (buffer_sizes[2] > 0)
        {
            if (extra_buffer_[7].count + buffer_sizes[2] > extra_buffer_[7].capacity)
            {
                extra_buffer_[7].resizeManaged(extra_buffer_[7].count + buffer_sizes[2] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[7].p+3*extra_buffer_[7].count, 3*buffer_sizes[2], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[7].q+3*extra_buffer_[7].count, 3*buffer_sizes[2], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[7].other+extra_buffer_[7].count, buffer_sizes[2], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            extra_buffer_[7].count += buffer_sizes[2];
        }
    }

    // second, odd ranks send downwards even ranks receive
    if (parallel.grid_rank()[1] % 2 == 1) // odd rank
    {
        buffer_sizes[0] = extra_buffer_[0].count;
        buffer_sizes[1] = extra_buffer_[3].count;
        buffer_sizes[2] = extra_buffer_[6].count;

        parallel.send_dim1(buffer_sizes, 3, parallel.grid_rank()[1]-1);
        if (buffer_sizes[0] > 0)
        {
            parallel.send_dim1(extra_buffer_[0].p, 3*extra_buffer_[0].count, parallel.grid_rank()[1]-1);
            parallel.send_dim1(extra_buffer_[0].q, 3*extra_buffer_[0].count, parallel.grid_rank()[1]-1);
            parallel.send_dim1(extra_buffer_[0].other, extra_buffer_[0].count, parallel.grid_rank()[1]-1);
        }
        if (buffer_sizes[1] > 0)
        {
            parallel.send_dim1(extra_buffer_[3].p, 3*extra_buffer_[3].count, parallel.grid_rank()[1]-1);
            parallel.send_dim1(extra_buffer_[3].q, 3*extra_buffer_[3].count, parallel.grid_rank()[1]-1);
            parallel.send_dim1(extra_buffer_[3].other, extra_buffer_[3].count, parallel.grid_rank()[1]-1);
        }
        if (buffer_sizes[2] > 0)
        {
            parallel.send_dim1(extra_buffer_[6].p, 3*extra_buffer_[6].count, parallel.grid_rank()[1]-1);
            parallel.send_dim1(extra_buffer_[6].q, 3*extra_buffer_[6].count, parallel.grid_rank()[1]-1);
            parallel.send_dim1(extra_buffer_[6].other, extra_buffer_[6].count, parallel.grid_rank()[1]-1);
        }
    }
    else
    {
        parallel.receive_dim1(buffer_sizes, 3, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);

        if (buffer_sizes[0] > 0)
        {
            if (extra_buffer_[1].count + buffer_sizes[0] > extra_buffer_[1].capacity)
            {
                extra_buffer_[1].resizeManaged(extra_buffer_[1].count + buffer_sizes[0] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[1].p+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].q+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].other+extra_buffer_[1].count, buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            extra_buffer_[1].count += buffer_sizes[0];
        }
        if (buffer_sizes[1] > 0)
        {
            if (extra_buffer_[4].count + buffer_sizes[1] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[1] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[4].p+3*extra_buffer_[4].count, 3*buffer_sizes[1], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[4].q+3*extra_buffer_[4].count, 3*buffer_sizes[1], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[4].other+extra_buffer_[4].count, buffer_sizes[1], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            extra_buffer_[4].count += buffer_sizes[1];
        }
        if (buffer_sizes[2] > 0)
        {
            if (extra_buffer_[7].count + buffer_sizes[2] > extra_buffer_[7].capacity)
            {
                extra_buffer_[7].resizeManaged(extra_buffer_[7].count + buffer_sizes[2] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[7].p+3*extra_buffer_[7].count, 3*buffer_sizes[2], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[7].q+3*extra_buffer_[7].count, 3*buffer_sizes[2], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[7].other+extra_buffer_[7].count, buffer_sizes[2], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            extra_buffer_[7].count += buffer_sizes[2];
        }
    }

    // third, even ranks send downwards odd ranks receive
    if (parallel.grid_rank()[1] % 2 == 0) // even rank
    {
        buffer_sizes[0] = extra_buffer_[0].count;
        buffer_sizes[1] = extra_buffer_[3].count;
        buffer_sizes[2] = extra_buffer_[6].count;

        parallel.send_dim1(buffer_sizes, 3, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
        if (buffer_sizes[0] > 0)
        {
            parallel.send_dim1(extra_buffer_[0].p, 3*extra_buffer_[0].count, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[0].q, 3*extra_buffer_[0].count, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[0].other, extra_buffer_[0].count, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
        }
        if (buffer_sizes[1] > 0)
        {
            parallel.send_dim1(extra_buffer_[3].p, 3*extra_buffer_[3].count, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[3].q, 3*extra_buffer_[3].count, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[3].other, extra_buffer_[3].count, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
        }
        if (buffer_sizes[2] > 0)
        {
            parallel.send_dim1(extra_buffer_[6].p, 3*extra_buffer_[6].count, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[6].q, 3*extra_buffer_[6].count, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[6].other, extra_buffer_[6].count, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
        }
    }
    else
    {
        parallel.receive_dim1(buffer_sizes, 3, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);

        if (buffer_sizes[0] > 0)
        {
            if (extra_buffer_[1].count + buffer_sizes[0] > extra_buffer_[1].capacity)
            {
                extra_buffer_[1].resizeManaged(extra_buffer_[1].count + buffer_sizes[0] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[1].p+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].q+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].other+extra_buffer_[1].count, buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            extra_buffer_[1].count += buffer_sizes[0];
        }
        if (buffer_sizes[1] > 0)
        {
            if (extra_buffer_[4].count + buffer_sizes[1] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[1] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[4].p+3*extra_buffer_[4].count, 3*buffer_sizes[1], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[4].q+3*extra_buffer_[4].count, 3*buffer_sizes[1], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[4].other+extra_buffer_[4].count, buffer_sizes[1], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            extra_buffer_[4].count += buffer_sizes[1];
        }
        if (buffer_sizes[2] > 0)
        {
            if (extra_buffer_[7].count + buffer_sizes[2] > extra_buffer_[7].capacity)
            {
                extra_buffer_[7].resizeManaged(extra_buffer_[7].count + buffer_sizes[2] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[7].p+3*extra_buffer_[7].count, 3*buffer_sizes[2], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[7].q+3*extra_buffer_[7].count, 3*buffer_sizes[2], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[7].other+extra_buffer_[7].count, buffer_sizes[2], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            extra_buffer_[7].count += buffer_sizes[2];
        }
    }

    // fourth, odd ranks send upwards even ranks receive
    if (parallel.grid_rank()[1] % 2 == 1) // odd rank
    {
        buffer_sizes[0] = extra_buffer_[2].count;
        buffer_sizes[1] = extra_buffer_[5].count;
        buffer_sizes[2] = extra_buffer_[8].count;

        parallel.send_dim1(buffer_sizes, 3, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
        if (buffer_sizes[0] > 0)
        {
            parallel.send_dim1(extra_buffer_[2].p, 3*extra_buffer_[2].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[2].q, 3*extra_buffer_[2].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[2].other, extra_buffer_[2].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
        }
        if (buffer_sizes[1] > 0)
        {
            parallel.send_dim1(extra_buffer_[5].p, 3*extra_buffer_[5].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[5].q, 3*extra_buffer_[5].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[5].other, extra_buffer_[5].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
        }
        if (buffer_sizes[2] > 0)
        {
            parallel.send_dim1(extra_buffer_[8].p, 3*extra_buffer_[8].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[8].q, 3*extra_buffer_[8].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.send_dim1(extra_buffer_[8].other, extra_buffer_[8].count, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
        }
    }
    else if (parallel.grid_rank()[1] > 0 || parallel.grid_size()[1] % 2 == 0)
    {
        parallel.receive_dim1(buffer_sizes, 3, (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);

        if (buffer_sizes[0] > 0)
        {
            if (extra_buffer_[1].count + buffer_sizes[0] > extra_buffer_[1].capacity)
            {
                extra_buffer_[1].resizeManaged(extra_buffer_[1].count + buffer_sizes[0] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[1].p+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].q+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].other+extra_buffer_[1].count, buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            extra_buffer_[1].count += buffer_sizes[0];
        }
        if (buffer_sizes[1] > 0)
        {
            if (extra_buffer_[4].count + buffer_sizes[1] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[1] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[4].p+3*extra_buffer_[4].count, 3*buffer_sizes[1], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[4].q+3*extra_buffer_[4].count, 3*buffer_sizes[1], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[4].other+extra_buffer_[4].count, buffer_sizes[1], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            extra_buffer_[4].count += buffer_sizes[1];
        }
        if (buffer_sizes[2] > 0)
        {
            if (extra_buffer_[7].count + buffer_sizes[2] > extra_buffer_[7].capacity)
            {
                extra_buffer_[7].resizeManaged(extra_buffer_[7].count + buffer_sizes[2] + extra_capacity_);
            }
            
            parallel.receive_dim1(extra_buffer_[7].p+3*extra_buffer_[7].count, 3*buffer_sizes[2], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[7].q+3*extra_buffer_[7].count, 3*buffer_sizes[2], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[7].other+extra_buffer_[7].count, buffer_sizes[2], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            extra_buffer_[7].count += buffer_sizes[2];
        }
    }

    // flush send buffer 4
    flushExtraBuffer(4);

    // now the remaining dimension needs to communicate
    // first, even ranks send upwards odd ranks receive
    if (parallel.grid_rank()[0] % 2 == 0) // even rank
    {
        buffer_sizes[0] = extra_buffer_[7].count;

        parallel.send_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            parallel.send_dim0(extra_buffer_[7].p, 3*extra_buffer_[7].count, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.send_dim0(extra_buffer_[7].q, 3*extra_buffer_[7].count, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.send_dim0(extra_buffer_[7].other, extra_buffer_[7].count, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
        }
    }
    
    if (parallel.grid_rank()[0] % 2 == 1 || (parallel.grid_rank()[0] == 0 && parallel.grid_size()[0] % 2 == 1)) // odd rank or rank 0 with odd grid size
    {
        parallel.receive_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            if (extra_buffer_[4].count + buffer_sizes[0] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[0] + extra_capacity_);
            }
            
            parallel.receive_dim0(extra_buffer_[4].p+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].q+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].other+extra_buffer_[4].count, buffer_sizes[0], (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            extra_buffer_[4].count += buffer_sizes[0];
        }
    }

    // second, odd ranks send downwards even ranks receive
    if (parallel.grid_rank()[0] % 2 == 1) // odd rank
    {
        buffer_sizes[0] = extra_buffer_[1].count;

        parallel.send_dim0(buffer_sizes, 1, parallel.grid_rank()[0]-1);

        if (buffer_sizes[0] > 0)
        {
            parallel.send_dim0(extra_buffer_[1].p, 3*extra_buffer_[1].count, parallel.grid_rank()[0]-1);
            parallel.send_dim0(extra_buffer_[1].q, 3*extra_buffer_[1].count, parallel.grid_rank()[0]-1);
            parallel.send_dim0(extra_buffer_[1].other, extra_buffer_[1].count, parallel.grid_rank()[0]-1);
        }
    }
    else
    {
        parallel.receive_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            if (extra_buffer_[4].count + buffer_sizes[0] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[0] + extra_capacity_);
            }
            
            parallel.receive_dim0(extra_buffer_[4].p+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].q+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].other+extra_buffer_[4].count, buffer_sizes[0], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            extra_buffer_[4].count += buffer_sizes[0];
        }
    }

    // third, even ranks send downwards odd ranks receive
    if (parallel.grid_rank()[0] % 2 == 0) // even rank
    {
        buffer_sizes[0] = extra_buffer_[1].count;

        parallel.send_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            parallel.send_dim0(extra_buffer_[1].p, 3*extra_buffer_[1].count, (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            parallel.send_dim0(extra_buffer_[1].q, 3*extra_buffer_[1].count, (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            parallel.send_dim0(extra_buffer_[1].other, extra_buffer_[1].count, (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
        }
    }
    else
    {
        parallel.receive_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            if (extra_buffer_[4].count + buffer_sizes[0] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[0] + extra_capacity_);
            }
            
            parallel.receive_dim0(extra_buffer_[4].p+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].q+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].other+extra_buffer_[4].count, buffer_sizes[0], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            extra_buffer_[4].count += buffer_sizes[0];
        }
    }

    // fourth, odd ranks send upwards even ranks receive
    if (parallel.grid_rank()[0] % 2 == 1) // odd rank
    {
        buffer_sizes[0] = extra_buffer_[7].count;

        parallel.send_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            parallel.send_dim0(extra_buffer_[7].p, 3*extra_buffer_[7].count, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.send_dim0(extra_buffer_[7].q, 3*extra_buffer_[7].count, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.send_dim0(extra_buffer_[7].other, extra_buffer_[7].count, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
        }
    }
    else if (parallel.grid_rank()[0] > 0 || parallel.grid_size()[0] % 2 == 0)
    {
        parallel.receive_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            if (extra_buffer_[4].count + buffer_sizes[0] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[0] + extra_capacity_);
            }
            
            parallel.receive_dim0(extra_buffer_[4].p+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].q+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].other+extra_buffer_[4].count, buffer_sizes[0], (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            extra_buffer_[4].count += buffer_sizes[0];
        }
    }

    nvtxRangePop();

    // ingest particles from send_buffer[4]

    nvtxRangePushA("moveParticles: ingest received particles");

    flushExtraBuffer(4);

    // periodic boundary conditions

    nvtxRangePushA("moveParticles: periodic boundary conditions");

    tripleReal * pos = (tripleReal *) p;

    thrust::for_each(thrust::device, pos, pos+num_particles_, [this] __host__ __device__ (tripleReal& a) {
        Real * pos_a = (Real *) &a;

        for (int i = 0; i < 3; i++)
        {
            if (pos_a[i] < 0.)
                pos_a[i] += boxSize_[i];
            if (pos_a[i] >= boxSize_[i])
                pos_a[i] -= boxSize_[i];
        }
    });

    nvtxRangePop();

#ifdef DEBUG_RADIX_SORT
    nvtxRangePushA("moveParticles: update row buffers");

    updateRowBuffers();

    nvtxRangePop();
#else
    nvtxRangePushA("moveParticles: compute rows for new particles");

    if (num_particles_ > start_capacity)
    {
        uint32_t * d_keys_temp = d_keys;

        success = cudaMallocAsync((void **) &d_keys, num_particles_ * sizeof(uint32_t), pcl_stream);
        if (success != cudaSuccess)
        {
            throw std::runtime_error("CUDA malloc failed in moveParticles");
        }

        success = cudaMemcpyAsync(d_keys, d_keys_temp, start_num_particles * sizeof(uint32_t), cudaMemcpyDeviceToDevice, pcl_stream);
        if (success != cudaSuccess)
        {
            throw std::runtime_error("CUDA memcpy failed in moveParticles");
        }

        cudaFreeAsync(d_keys_temp, pcl_stream);
    }

    if (num_particles_ > start_num_particles)
        compute_rows<<<(num_particles_-start_num_particles+127)/128, 128, 0, pcl_stream>>>(this, d_keys, start_num_particles);

    nvtxRangePop();

    nvtxRangePushA("moveParticles: sort particles");

    unsigned long long * d_indices_in;
    success = cudaMallocAsync((void **) &d_indices_in, num_particles_ * sizeof(unsigned long long) * 2L, pcl_stream);
    if (success != cudaSuccess)
    {
        throw std::runtime_error("CUDA malloc failed in moveParticles");
    }

    computeSortIndices(d_keys, d_indices_in, &d_temp, pcl_stream);

    nvtxRangePop();

    reorderParticles(d_indices_in, d_temp, pcl_stream);

    cudaFreeAsync(d_temp, pcl_stream);
    cudaFreeAsync(d_indices_in, pcl_stream);

    update_pointers<<<(num_row_buffers_+137)/128, 128, 0, pcl_stream>>>(this, nullptr);

    success = cudaStreamSynchronize(pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA kernel reorder_data");
    }

    cudaStreamDestroy(pcl_stream);
#endif

    nvtxRangePop(); // moveParticles
}

// project particles
template <typename part, typename part_info>
__global__ void project_particles(perfParticles<part, part_info> * pcl, Real * target, int projection_order, long start_idx, const long * jump, int stencil_k)
{
    int row = blockIdx.x;
    int thread_id = threadIdx.x;

    for (int idx = thread_id; idx < pcl->row_buffers_[row].count; idx += 128)
    {
        pcl->project_particle(target, projection_order, start_idx, jump, row, stencil_k, idx);
    }
}

// density projection
template <typename part, typename part_info>
void perfParticles<part, part_info>::meshprojection_project(Field<Real> * target, int projection_order)
{
    long start_idx = target->lattice().siteFirst();
    long jump[3];

    Real * data = target->data();

    for (int i = 0; i < 3; i++)
    {
        jump[i] = target->lattice().jump(i);
    }

    // Loop over the projection stencil
    for (int k = 0; k < 4; k++)
    {
        project_particles<<<num_row_buffers_, 128>>>(this, data, projection_order, start_idx, jump, k);
        auto status = cudaDeviceSynchronize();

        if (status != cudaSuccess)
        {
            std::cerr << "CUDA kernel failed: " << cudaGetErrorString(status) << std::endl;
            throw std::runtime_error("Error in CUDA kernel project_particles");
        }
    }
}

// project particle
template <typename part, typename part_info>
__device__ void perfParticles<part, part_info>::project_particle(Real * target, int projection_order, long start_idx, const long * jump, int row, int stencil_k, int idx)
{
    int coord[3];
    Real frac[3];
    Real dx = boxSize_[0] / lat_size_[0];

    Real weight = part_global_info_.mass / (dx * dx * dx);

    coord[1] = static_cast<int>(floor(row_buffers_[row].p[1]*lat_size_[1]));
    coord[2] = static_cast<int>(floor(row_buffers_[row].p[2]*lat_size_[2]));
    coord[0] = 0;

    long site = start_idx + coord[0]*jump[0] + (coord[1] - coordSkip_[1])*jump[1] + (coord[2] - coordSkip_[0])*jump[2];

    if (stencil_k % 2)
    {
        site += jump[1];
        coord[1]++;
    }
    if (stencil_k / 2)
    {
        site += jump[2];
        coord[2]++;
    }

    Real * target_row = target+site;

    if (projection_order == PERFPARTICLES_CIC)
    {
        coord[0] = static_cast<int>(floor(row_buffers_[row].p[3*idx]*lat_size_[0]));
        
        frac[0] = Real(1)-(row_buffers_[row].p[3*idx]*lat_size_[0] - coord[0]);
        frac[1] = Real(1)-fabs(row_buffers_[row].p[3*idx+1]*lat_size_[1] - coord[1]);
        frac[2] = Real(1)-fabs(row_buffers_[row].p[3*idx+2]*lat_size_[2] - coord[2]);

        atomicAdd(target_row+coord[0], weight * frac[0] * frac[1] * frac[2]);
        atomicAdd(target_row+(coord[0]+1) % lat_size_[0], weight * (Real(1)-frac[0]) * frac[1] * frac[2]);
    }
    else if (projection_order == PERFPARTICLES_NGP)
    {
        coord[0] = (int) floor(row_buffers_[row].p[3*idx]*lat_size_[0]);
        
        frac[0] = (row_buffers_[row].p[3*idx]*lat_size_[0] - coord[0]) < Real(0.5) ? Real(1) : Real(0);
        frac[1] = fabs(row_buffers_[row].p[3*idx+1]*lat_size_[1] - coord[1]) < Real(0.5) ? Real(1) : Real(0);
        frac[2] = fabs(row_buffers_[row].p[3*idx+2]*lat_size_[2] - coord[2]) < Real(0.5) ? Real(1) : Real(0);

        target_row[coord[0]] += weight * frac[0] * frac[1] * frac[2];
        target_row[(coord[0]+1) % lat_size_[0]] += weight * (Real(1)-frac[0]) * frac[1] * frac[2];
    }
}

#endif
