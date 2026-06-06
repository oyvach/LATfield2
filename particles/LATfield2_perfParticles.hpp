#ifndef LATFIELD2_PERFPARTICLES_HPP
#define LATFIELD2_PERFPARTICLES_HPP

#include <algorithm>
#include <cstdlib>
#include <cstring>

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
    void resize(int new_capacity, bool managed = false);
    void resizeManaged(int new_capacity);
};

// resize function
template <typename p_type, typename q_type, typename other_type>
void row_buffer<p_type, q_type, other_type>::resize(int new_capacity, bool managed)
{
    if (capacity < 0)
    {
        throw std::runtime_error("Trying to resize unmanaged row buffer");
    }
    p_type * new_p;
    q_type * new_q;
    other_type * new_other;
    // Reallocate the arrays using realloc
    if (!managed)
    {
        new_p = (p_type *)realloc(p, new_capacity * 3 * sizeof(p_type));
        new_q = (q_type *)realloc(q, new_capacity * 3 * sizeof(q_type));
        new_other = (other_type *)realloc(other, new_capacity * sizeof(other_type));
    }
    else
    {
        new_p = cuda_reallocManaged(p, capacity * 3, new_capacity * 3);
        new_q = cuda_reallocManaged(q, capacity * 3, new_capacity * 3);
        new_other = cuda_reallocManaged(other, capacity, new_capacity);
    }
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
    resize(new_capacity, true);
}

template <typename part, typename part_info>
class perfParticles;

template <typename part, typename part_info>
__global__ void update_pointers(perfParticles<part, part_info> * pcl, unsigned long long int * send_begin);

template <typename part, typename part_info>
__global__ void compute_rows(perfParticles<part, part_info> * pcl, uint32_t * row, unsigned long long int starting_idx = 0);

__global__ void compute_row_offsets(const uint32_t * keys, unsigned long long int num_particles, unsigned long long int * offsets, uint32_t num_rows);

template <typename T>
__global__ void reorder_data(unsigned long long int * indices_out, T * data_in, T * data_out, size_t stride, size_t ndata);

template <typename part, typename part_info, typename UpdateFunct>
__global__ void update_particles(perfParticles<part, part_info> * pcl, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, Real * v2, bool copyout = true, void ** vparams = nullptr);

template <typename part, typename part_info>
__global__ void project_particles(perfParticles<part, part_info> * pcl, Real * target, int projection_order, long start_idx, long jump0, long jump1, long jump2, int stencil_k);

// Particle handler class
template <typename part, typename part_info>
class perfParticles
{
    public:
        // Constructor
        perfParticles(bool use_managed = false);
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
        size_t temporaryWorkspaceBytes() const;
        static size_t temporaryWorkspaceBytesForParticles(uint64_t particle_capacity, uint32_t num_row_buffers);
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
        Real updateVel(UpdateFunct updateVel_funct, double dtau, Field<Real> ** fields = NULL, int nfields = 0, double * params = NULL, double * output = NULL, int * reduce_type = NULL, int noutput = 0, void **vparams = NULL);
        // Move particles
        template <typename UpdateFunct>
        void moveParticles(UpdateFunct move_funct, double dtau, Field<Real> ** fields = NULL, int nfields = 0, double * params = NULL, double * output = NULL, int * reduce_type = NULL, int noutput = 0, void **vparams = NULL);
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

        bool managed_runtime_;

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

            if (managed_runtime_)
            {
                p = cuda_reallocManaged(p, total_capacity_ * 3, new_total_capacity * 3);
                q = cuda_reallocManaged(q, total_capacity_ * 3, new_total_capacity * 3);
                other = cuda_reallocManaged(other, total_capacity_, new_total_capacity);
            }
            else
            {
                p = cuda_realloc(p, total_capacity_ * 3, new_total_capacity * 3);
                q = cuda_realloc(q, total_capacity_ * 3, new_total_capacity * 3);
                other = cuda_realloc(other, total_capacity_, new_total_capacity);
            }

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
        void prepareComm(unsigned long long int * send_begin, uint32_t ** d_keys, void ** d_temp, bool * d_temp_private, cudaStream_t & stream);

        // helper function for sorting particles
        void computeSortIndices(uint32_t * d_keys, unsigned long long int * d_indices, void ** d_temp, bool * d_temp_private, cudaStream_t & stream);

        // helper function to reorder particles
        void reorderParticles(unsigned long long int * d_indices, void * d_temp, cudaStream_t & stream);

        void acquireTemporaryWorkspace(void ** d_temp, bool * d_temp_private, size_t bytes, cudaStream_t & stream, const char * context);
        void releaseTemporaryWorkspace(void ** d_temp, bool * d_temp_private, cudaStream_t & stream);

        // helper function for particle-mesh projection
        __device__ void project_particle(Real * target, int projection_order, long start_idx, long jump0, long jump1, long jump2, int row, int stencil_k, int idx);

        // helper function for particle updates
        template <typename UpdateFunct>
        void updateParticles(UpdateFunct update_funct, double dtau, Field<Real> ** fields = NULL, int nfields = 0, double * params = NULL, double * output = NULL, int * reduce_type = NULL, int noutput = 0, Real * maxvel = nullptr, bool copyout = true, bool async = false, void ** vparams = NULL);

        template <typename UpdateFunct>
        __host__ __device__ auto updateParticle(int row, int idx, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int noutput, bool copyout = true, void ** vparams = NULL);

        // friend functions

        template <typename part2, typename part_info2>
        friend __global__ void update_pointers(perfParticles<part2, part_info2> * pcl, unsigned long long int * send_begin);

        template <typename part2, typename part_info2>
        friend __global__ void compute_rows(perfParticles<part2, part_info2> * pcl, uint32_t * row, unsigned long long int starting_idx);

        template <typename part2, typename part_info2, typename UpdateFunct>
        friend __global__ void update_particles(perfParticles<part2, part_info2> * pcl, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, Real * v2, bool copyout, void ** vparams);

        template <typename part2, typename part_info2>
        friend __global__ void project_particles(perfParticles<part2, part_info2> * pcl, Real * target, int projection_order, long start_idx, long jump0, long jump1, long jump2, int stencil_k);
};


// Constructor
template <typename part, typename part_info>
perfParticles<part, part_info>::perfParticles(bool use_managed)
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
    managed_runtime_ = use_managed;
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
        if (managed_runtime_)
        {
            cudaFree(extra_buffer_[i].p);
            cudaFree(extra_buffer_[i].q);
            cudaFree(extra_buffer_[i].other);
        }
        else
        {
            free(extra_buffer_[i].p);
            free(extra_buffer_[i].q);
            free(extra_buffer_[i].other);
        }
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

    if (managed_runtime_)
    {
        success = cudaMallocManaged(&p, total_capacity_ * 3 * sizeof(Real));
    }
    else
    {
        success = cudaMalloc(&p, total_capacity_ * 3 * sizeof(Real)); // this should be on device
    }

    if (success != cudaSuccess)
    {
        std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Memory allocation failed for global position array in perfParticles::initialize");
    }

    if (managed_runtime_)
    {
        success = cudaMallocManaged(&q, total_capacity_ * 3 * sizeof(Real));
    }
    else
    {
        success = cudaMalloc(&q, total_capacity_ * 3 * sizeof(Real)); // this should be on device
    }

    if (success != cudaSuccess)
    {
        std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Memory allocation failed for global momentum array in perfParticles::initialize");
    }

    if (managed_runtime_)
    {
        success = cudaMallocManaged(&other, total_capacity_ * sizeof(long));
    }
    else
    {
        success = cudaMalloc(&other, total_capacity_ * sizeof(long)); // this should be on device
    }

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
template <typename T>
__global__ void reorder_data(unsigned long long int * indices_out, T * data_in, T * data_out, size_t stride, size_t ndata)
{
    unsigned long long int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < ndata)
    {
        unsigned long long int pcl_idx = indices_out[idx / stride];

        data_out[idx] = data_in[pcl_idx * stride + (idx % stride)];
    }
}

template <typename part, typename part_info>
size_t perfParticles<part, part_info>::temporaryWorkspaceBytes() const
{
    return temporaryWorkspaceBytesForParticles(total_capacity_, num_row_buffers_);
}

template <typename part, typename part_info>
size_t perfParticles<part, part_info>::temporaryWorkspaceBytesForParticles(uint64_t particle_capacity, uint32_t num_row_buffers)
{
    if (particle_capacity == 0) return 0;

    size_t workspace_bytes = particle_capacity * 3 * sizeof(Real);
    size_t temp_storage_bytes = 0;
    uint32_t * keys32 = nullptr;
    unsigned long long int * indices = nullptr;

    int comm_end_bit = static_cast<int>(ceil(log2(static_cast<double>(num_row_buffers + 10))));
    cub::DeviceRadixSort::SortPairs(nullptr, temp_storage_bytes, keys32, keys32, indices, indices, particle_capacity, 0, comm_end_bit);
    workspace_bytes = std::max(workspace_bytes, temp_storage_bytes);

    temp_storage_bytes = 0;
    int row_end_bit = (num_row_buffers > 1) ? static_cast<int>(ceil(log2(static_cast<double>(num_row_buffers)))) : 0;
    cub::DeviceRadixSort::SortPairs(nullptr, temp_storage_bytes, keys32, keys32, indices, indices, particle_capacity, 0, row_end_bit);
    workspace_bytes = std::max(workspace_bytes, temp_storage_bytes);

    temp_storage_bytes = 0;
    unsigned long long int * row_offsets = nullptr;
#ifdef SINGLE
    cub::DeviceSegmentedSort::SortPairs(nullptr, temp_storage_bytes, keys32, keys32, indices, indices, particle_capacity, num_row_buffers, row_offsets, row_offsets);
#else
    uint64_t * keys64 = nullptr;
    cub::DeviceSegmentedSort::SortPairs(nullptr, temp_storage_bytes, keys64, keys64, indices, indices, particle_capacity, num_row_buffers, row_offsets, row_offsets);
#endif
    workspace_bytes = std::max(workspace_bytes, temp_storage_bytes);

    return workspace_bytes;
}

template <typename part, typename part_info>
void perfParticles<part, part_info>::acquireTemporaryWorkspace(void ** d_temp, bool * d_temp_private, size_t bytes, cudaStream_t & stream, const char * context)
{
    if (bytes == 0)
    {
        releaseTemporaryWorkspace(d_temp, d_temp_private, stream);
        return;
    }

    releaseTemporaryWorkspace(d_temp, d_temp_private, stream);

#ifdef FFT3D
    nvtxRangePushA("perfParticles: acquire shared temporary workspace");
    tempMemory.reserveDeviceWorkspaceBytes(bytes, context);
    *d_temp = tempMemory.deviceWorkspace();
    *d_temp_private = false;
    nvtxRangePop();
#else
    nvtxRangePushA("perfParticles: allocate private temporary workspace");
    auto success = cudaMallocAsync(d_temp, bytes, stream);
    if (success != cudaSuccess)
    {
        nvtxRangePop();
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for particle temporary workspace");
    }
    *d_temp_private = true;
    nvtxRangePop();
#endif
}

template <typename part, typename part_info>
void perfParticles<part, part_info>::releaseTemporaryWorkspace(void ** d_temp, bool * d_temp_private, cudaStream_t & stream)
{
    if (*d_temp != nullptr && *d_temp_private)
    {
        cudaFreeAsync(*d_temp, stream);
    }
    *d_temp = nullptr;
    *d_temp_private = false;
}

// prepare communication
template <typename part, typename part_info>
void perfParticles<part, part_info>::prepareComm(unsigned long long int * send_begin, uint32_t ** d_keys, void ** d_temp, bool * d_temp_private, cudaStream_t & stream)
{
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

    acquireTemporaryWorkspace(d_temp, d_temp_private, temp_storage_bytes, stream, "perfParticles::prepareComm");

    // sort
    cub::DeviceRadixSort::SortPairs(*d_temp, temp_storage_bytes, d_keys_in, *d_keys, d_indices_in, d_indices_out, num_particles_, 0, end_bit, stream);
    nvtxRangePop();

    cudaFreeAsync(d_keys_in, stream);

    reorderParticles(d_indices_out, *d_temp, stream);

    cudaFreeAsync(d_indices_in, stream);

    nvtxRangePushA("prepareComm: update pointers");

    update_pointers<<<(num_row_buffers_+137)/128, 128, 0, stream>>>(this, send_begin);

    /*success = cudaStreamSynchronize(stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA kernel update_pointers");
    }*/
    nvtxRangePop();

    rows_sorted_ = false;
}

// compute sort indices
template <typename part, typename part_info>
void perfParticles<part, part_info>::computeSortIndices(uint32_t * d_keys, unsigned long long int * d_indices, void ** d_temp, bool * d_temp_private, cudaStream_t & stream)
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

    unsigned long long int * row_offsets_query = nullptr;
#ifdef SINGLE
    uint32_t * xkeys_query = nullptr;
#else
    uint64_t * xkeys_query = nullptr;
#endif
    cub::DeviceSegmentedSort::SortPairs(nullptr, temp_storage_bytes2, xkeys_query, xkeys_query, d_indices_out, d_indices, num_particles_, num_row_buffers_, row_offsets_query, row_offsets_query, stream);
    temp_storage_bytes = std::max(temp_storage_bytes, temp_storage_bytes2);

    acquireTemporaryWorkspace(d_temp, d_temp_private, temp_storage_bytes, stream, "perfParticles::computeSortIndices");

    // sort by row key
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
    reorder_data<<<(num_particles_ + 127)/128, 128, 0, stream>>>(d_indices, other, (long *) d_temp, 1, num_particles_);

    // copy back to other
    auto success = cudaMemcpyAsync(other, d_temp, num_particles_ * sizeof(long), cudaMemcpyDeviceToDevice, stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA memcpy failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA memcpy for other array");
    }

    // reorder position data
    reorder_data<<<(num_particles_ * 3 + 127)/128, 128, 0, stream>>>(d_indices, p, (Real*) d_temp, 3, num_particles_ * 3);

    // reorder momentum data using old position array as target
    reorder_data<<<(num_particles_ * 3 + 127)/128, 128, 0, stream>>>(d_indices, q, p, 3, num_particles_ * 3);

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

// update row buffers
template <typename part, typename part_info>
void perfParticles<part, part_info>::updateRowBuffers(unsigned long long int * send_begin)
{
#ifndef DEBUG_RADIX_SORT
    uint32_t * d_keys_in = nullptr;
    unsigned long long int * d_indices = nullptr;
    void * d_temp = nullptr;
    bool d_temp_private = false;
    cudaStream_t pcl_stream;

    auto success = cudaStreamCreateWithFlags(&pcl_stream, cudaStreamNonBlocking);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA stream creation failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA stream creation for moveParticles");
    }
#endif

    // check if adding buffer is empty, otherwise add particles to global buffer
    nvtxRangePushA("updateRowBuffers: add particles");

    flushExtraBuffer(4);

    nvtxRangePop();

    nvtxRangePushA("updateRowBuffers: sort");
 #ifdef DEBUG_RADIX_SORT
    // tuple-pointer to particle properties
    tripleReal * pos = (tripleReal *) p;
    tripleReal * vel = (tripleReal *) q;

    // create zip iterator from tuples
    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(pos, vel, other));
    auto zip_end = thrust::make_zip_iterator(thrust::make_tuple(pos + num_particles_, vel + num_particles_, other + num_particles_));

    using zip_iter = typename decltype(zip_begin)::value_type;

    // sort particles
    thrust::sort(thrust::device, zip_begin, zip_end, [this] __device__ (zip_iter a, zip_iter b) {
        auto pos_a = thrust::get<0>(a);
        auto pos_b = thrust::get<0>(b);

        int row_a = computeRow(((Real *) (&pos_a))[1], ((Real *) (&pos_a))[2]);
        int row_b = computeRow(((Real *) (&pos_b))[1], ((Real *) (&pos_b))[2]);

        if (row_a == row_b)
        {
            return ((Real *) (&pos_a))[0] < ((Real *) (&pos_b))[0];
        }

        return row_a < row_b;
    });
#else
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

    computeSortIndices(d_keys_in, d_indices, &d_temp, &d_temp_private, pcl_stream);

    reorderParticles(d_indices, d_temp, pcl_stream);

    releaseTemporaryWorkspace(&d_temp, &d_temp_private, pcl_stream);
    cudaFreeAsync(d_indices, pcl_stream);

    success = cudaStreamSynchronize(pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error after cudaStreamSynchronize in updateRowBuffers");
    }

    cudaStreamDestroy(pcl_stream);
#endif
    nvtxRangePop();

    nvtxRangePushA("updateRowBuffers: update pointers");

    update_pointers<<<(num_row_buffers_+137)/128, 128>>>(this, send_begin);

    auto success2 = cudaDeviceSynchronize();

    if (success2 != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success2) << std::endl;
        throw std::runtime_error("Error after cudaDeviceSynchronize in updateRowBuffers");
    }

    nvtxRangePop();

    rows_sorted_ = true;
}

// update particles
template <typename part, typename part_info, typename UpdateFunct>
__global__ void update_particles(perfParticles<part, part_info> * pcl, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, Real * v2, bool copyout, void ** vparams)
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

    using return_type = decltype(pcl->updateParticle(row, thread_id, update_funct, dtau, fields, nfields, params, output_thread, noutput, copyout, vparams));

    for (int idx = thread_id; idx < pcl->row_buffers_[row].count; idx += 128)
    {
        if constexpr (std::is_same_v<return_type, void>)
        {
            pcl->updateParticle(row, idx, update_funct, dtau, fields, nfields, params, output_pcl, noutput, copyout, vparams);
        }
        else
        {
            Real v2_pcl = pcl->updateParticle(row, idx, update_funct, dtau, fields, nfields, params, output_pcl, noutput, copyout, vparams);

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
void perfParticles<part, part_info>::updateParticles(UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, Real * maxvel, bool copyout, bool async, void ** vparams)
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
    Real * d_v2 = nullptr;

    Field<Real> ** d_fields = nullptr;
    double * d_output = nullptr;
    int * d_reduce_type = nullptr;

    if (maxvel != nullptr)
    {
        auto successv2 = cudaMalloc(&d_v2, sizeof(Real));

        if (successv2 != cudaSuccess)
        {
            std::cerr << "cudaMalloc failed: " << cudaGetErrorString(successv2) << std::endl;
            throw std::runtime_error("Memory allocation failed for max velocity in perfParticles::updateParticles");
        }

        successv2 = cudaMemcpy(d_v2, &v2, sizeof(Real), cudaMemcpyDefault);

        if (successv2 != cudaSuccess)
        {
            std::cerr << "cudaMemcpy failed: " << cudaGetErrorString(successv2) << std::endl;
            throw std::runtime_error("Memory copy failed for max velocity in perfParticles::updateParticles");
        }
    }

    if (nfields > 0)
    {
        auto successfields = cudaMalloc(&d_fields, nfields * sizeof(Field<Real> *));

        if (successfields != cudaSuccess)
        {
            std::cerr << "cudaMalloc failed: " << cudaGetErrorString(successfields) << std::endl;
            throw std::runtime_error("Memory allocation failed for field pointers in perfParticles::updateParticles");
        }

        successfields = cudaMemcpy(d_fields, fields, nfields * sizeof(Field<Real> *), cudaMemcpyDefault);

        if (successfields != cudaSuccess)
        {
            std::cerr << "cudaMemcpy failed: " << cudaGetErrorString(successfields) << std::endl;
            throw std::runtime_error("Memory copy failed for field pointers in perfParticles::updateParticles");
        }
    }

    if (noutput > 0)
    {
        auto successoutput = cudaMalloc(&d_output, noutput * sizeof(double));

        if (successoutput != cudaSuccess)
        {
            std::cerr << "cudaMalloc failed: " << cudaGetErrorString(successoutput) << std::endl;
            throw std::runtime_error("Memory allocation failed for output array in perfParticles::updateParticles");
        }

        successoutput = cudaMemcpy(d_output, output, noutput * sizeof(double), cudaMemcpyDefault);

        if (successoutput != cudaSuccess)
        {
            std::cerr << "cudaMemcpy failed: " << cudaGetErrorString(successoutput) << std::endl;
            throw std::runtime_error("Memory copy failed for output array in perfParticles::updateParticles");
        }

        successoutput = cudaMalloc(&d_reduce_type, noutput * sizeof(int));

        if (successoutput != cudaSuccess)
        {
            std::cerr << "cudaMalloc failed: " << cudaGetErrorString(successoutput) << std::endl;
            throw std::runtime_error("Memory allocation failed for reduce type array in perfParticles::updateParticles");
        }

        successoutput = cudaMemcpy(d_reduce_type, reduce_type, noutput * sizeof(int), cudaMemcpyDefault);

        if (successoutput != cudaSuccess)
        {
            std::cerr << "cudaMemcpy failed: " << cudaGetErrorString(successoutput) << std::endl;
            throw std::runtime_error("Memory copy failed for reduce type array in perfParticles::updateParticles");
        }
    }

    update_particles<<<num_row_buffers_, 128>>>(this, update_funct, dtau, d_fields, nfields, params, d_output, d_reduce_type, noutput, d_v2, copyout, vparams);

    if (nfields > 0)
    {
        cudaFreeAsync(d_fields, 0);
    }

    if (noutput > 0)
    {
        cudaMemcpyAsync(output, d_output, noutput * sizeof(double), cudaMemcpyDefault);
        cudaFreeAsync(d_output, 0);
        cudaFreeAsync(d_reduce_type, 0);
    }

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
        auto successcpy = cudaMemcpy(&v2, d_v2, sizeof(Real), cudaMemcpyDefault);
        if (successcpy != cudaSuccess)
        {
            std::cerr << "cudaMemcpy failed: " << cudaGetErrorString(successcpy) << std::endl;
            throw std::runtime_error("Memory copy failed for max velocity in perfParticles::updateParticles");
        }

        cudaFree(d_v2);

        *maxvel = sqrt(v2);
    }
}

// Particle update helper function
template <typename part, typename part_info>
template <typename UpdateFunct>
__host__ __device__ auto perfParticles<part, part_info>::updateParticle(int row, int idx, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int noutput, bool copyout, void ** vparams)
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

    using return_type = decltype(update_funct(dtau, dx, &pcl, frac, part_global_info_, fields, sites, nfields, params, output, noutput, vparams));

    if constexpr (std::is_same_v<return_type, void>)
    {
        update_funct(dtau, dx, &pcl, frac, part_global_info_, fields, sites, nfields, params, output, noutput, vparams);
    }
    else
    {
        v2 = update_funct(dtau, dx, &pcl, frac, part_global_info_, fields, sites, nfields, params, output, noutput, vparams);
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
Real perfParticles<part, part_info>::updateVel(UpdateFunct updateVel_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, void ** vparams)
{
    Real maxvel = 0.;

    updateParticles(updateVel_funct, dtau, fields, nfields, params, output, reduce_type, noutput, &maxvel, true, false, vparams);

    return maxvel;
}

// Move particles
template <typename part, typename part_info>
template <typename UpdateFunct>
void perfParticles<part, part_info>::moveParticles(UpdateFunct move_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, void ** vparams)
{
    nvtxRangePushA("moveParticles");

    unsigned long long int send_begin[10];
    unsigned long long int * d_send_begin;
#ifndef DEBUG_RADIX_SORT
    uint32_t * d_keys = nullptr;
    void * d_temp = nullptr;
    bool d_temp_private = false;
    uint64_t start_num_particles;
    uint64_t start_capacity = total_capacity_;
    cudaStream_t pcl_stream;

    auto success = cudaStreamCreateWithFlags(&pcl_stream, cudaStreamNonBlocking);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA stream creation failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA stream creation for moveParticles");
    }

    success = cudaMallocAsync(&d_send_begin, 10 * sizeof(unsigned long long int), pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_send_begin");
    }
#else
    auto success = cudaMalloc(&d_send_begin, 10 * sizeof(unsigned long long int));
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA malloc for d_send_begin");
    }
#endif

    nvtxRangePushA("moveParticles: updateParticles");

    updateParticles(move_funct, dtau, fields, nfields, params, output, reduce_type, noutput, nullptr, true, false, vparams);

    nvtxRangePop();

    thrust::for_each(thrust::device, row_buffers_, row_buffers_ + num_row_buffers_, [] __device__ (row_buffer<Real, Real, long> & rb) { rb.sorted = false;});

    rows_sorted_ = false;

#ifdef DEBUG_RADIX_SORT
    nvtxRangePushA("moveParticles: updateRowBuffers");

    updateRowBuffers(d_send_begin);

    cudaMemcpy(send_begin, d_send_begin, 10 * sizeof(unsigned long long int), cudaMemcpyDefault);
    cudaFree(d_send_begin);

    nvtxRangePop();
#else
    nvtxRangePushA("moveParticles: prepareComm");

    prepareComm(d_send_begin, &d_keys, &d_temp, &d_temp_private, pcl_stream);

    releaseTemporaryWorkspace(&d_temp, &d_temp_private, pcl_stream);
    cudaMemcpyAsync(send_begin, d_send_begin, 10 * sizeof(unsigned long long int), cudaMemcpyDefault, pcl_stream);
    cudaFreeAsync(d_send_begin, pcl_stream);

    success = cudaStreamSynchronize(pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error after cudaStreamSynchronize in moveParticles");
    }

    nvtxRangePop();
#endif
    // memory movement

    num_particles_ = send_begin[0];
#ifndef DEBUG_RADIX_SORT
    start_num_particles = num_particles_;
#endif

    auto env_var_enabled = [](const char* name)
    {
        const char* value = std::getenv(name);
        if (value == nullptr) return false;
        return !(strcmp(value, "0") == 0 || strcmp(value, "false") == 0 || strcmp(value, "FALSE") == 0);
    };

    const bool cuda_aware_mpi_hint =
        env_var_enabled("LATFIELD2_ENABLE_CUDA_AWARE_MPI") ||
        env_var_enabled("MPICH_GPU_SUPPORT_ENABLED") ||
        env_var_enabled("MV2_USE_CUDA") ||
        env_var_enabled("PSM2_CUDA") ||
        env_var_enabled("OMPI_MCA_opal_cuda_support") ||
        env_var_enabled("OMPI_MCA_mpi_cuda_support");
    const bool cuda_aware_mpi_active =
        !env_var_enabled("LATFIELD2_DISABLE_CUDA_AWARE_MPI") &&
        cuda_aware_mpi_hint;

    if (cuda_aware_mpi_active)
    {
        struct device_particle_buffer
        {
            Real * p;
            Real * q;
            long * other;
            uint64_t count;
            uint64_t capacity;
        };

        auto class_count = [&](int proc) -> long
        {
            return static_cast<long>(send_begin[proc + 1] - send_begin[proc]);
        };

        auto send_dim1_sizes = [&](int proc0, int proc1, int proc2, int to)
        {
            long sizes[3] = {class_count(proc0), class_count(proc1), class_count(proc2)};
            parallel.send_dim1(sizes, 3, to);
        };

        long dim1_recv_sizes[4][3] = {};

        nvtxRangePushA("moveParticles: CUDA-aware particle size preflight");

        if (parallel.grid_rank()[1] % 2 == 0)
        {
            send_dim1_sizes(2, 5, 8, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
        }
        if (parallel.grid_rank()[1] % 2 == 1 || (parallel.grid_rank()[1] == 0 && parallel.grid_size()[1] % 2 == 1))
        {
            parallel.receive_dim1(dim1_recv_sizes[0], 3, (parallel.grid_rank()[1] + parallel.grid_size()[1] - 1) % parallel.grid_size()[1]);
        }

        if (parallel.grid_rank()[1] % 2 == 1)
        {
            send_dim1_sizes(0, 3, 6, parallel.grid_rank()[1] - 1);
        }
        else
        {
            parallel.receive_dim1(dim1_recv_sizes[1], 3, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
        }

        if (parallel.grid_rank()[1] % 2 == 0)
        {
            send_dim1_sizes(0, 3, 6, (parallel.grid_rank()[1] + parallel.grid_size()[1] - 1) % parallel.grid_size()[1]);
        }
        else
        {
            parallel.receive_dim1(dim1_recv_sizes[2], 3, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
        }

        if (parallel.grid_rank()[1] % 2 == 1)
        {
            send_dim1_sizes(2, 5, 8, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
        }
        else if (parallel.grid_rank()[1] > 0 || parallel.grid_size()[1] % 2 == 0)
        {
            parallel.receive_dim1(dim1_recv_sizes[3], 3, (parallel.grid_rank()[1] + parallel.grid_size()[1] - 1) % parallel.grid_size()[1]);
        }

        uint64_t route1_capacity = static_cast<uint64_t>(class_count(1));
        uint64_t route7_capacity = static_cast<uint64_t>(class_count(7));
        uint64_t dim1_final_capacity = 0;
        for (int phase = 0; phase < 4; phase++)
        {
            route1_capacity += static_cast<uint64_t>(dim1_recv_sizes[phase][0]);
            dim1_final_capacity += static_cast<uint64_t>(dim1_recv_sizes[phase][1]);
            route7_capacity += static_cast<uint64_t>(dim1_recv_sizes[phase][2]);
        }

        auto send_dim0_size = [&](long size, int to)
        {
            parallel.send_dim0(&size, 1, to);
        };

        long dim0_recv_sizes[4] = {};

        if (parallel.grid_rank()[0] % 2 == 0)
        {
            send_dim0_size(static_cast<long>(route7_capacity), (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
        }
        if (parallel.grid_rank()[0] % 2 == 1 || (parallel.grid_rank()[0] == 0 && parallel.grid_size()[0] % 2 == 1))
        {
            parallel.receive_dim0(&dim0_recv_sizes[0], 1, (parallel.grid_rank()[0] + parallel.grid_size()[0] - 1) % parallel.grid_size()[0]);
        }

        if (parallel.grid_rank()[0] % 2 == 1)
        {
            send_dim0_size(static_cast<long>(route1_capacity), parallel.grid_rank()[0] - 1);
        }
        else
        {
            parallel.receive_dim0(&dim0_recv_sizes[1], 1, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
        }

        if (parallel.grid_rank()[0] % 2 == 0)
        {
            send_dim0_size(static_cast<long>(route1_capacity), (parallel.grid_rank()[0] + parallel.grid_size()[0] - 1) % parallel.grid_size()[0]);
        }
        else
        {
            parallel.receive_dim0(&dim0_recv_sizes[2], 1, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
        }

        if (parallel.grid_rank()[0] % 2 == 1)
        {
            send_dim0_size(static_cast<long>(route7_capacity), (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
        }
        else if (parallel.grid_rank()[0] > 0 || parallel.grid_size()[0] % 2 == 0)
        {
            parallel.receive_dim0(&dim0_recv_sizes[3], 1, (parallel.grid_rank()[0] + parallel.grid_size()[0] - 1) % parallel.grid_size()[0]);
        }

        uint64_t dim0_final_capacity = 0;
        for (int phase = 0; phase < 4; phase++)
        {
            dim0_final_capacity += static_cast<uint64_t>(dim0_recv_sizes[phase]);
        }

        nvtxRangePop();

        uint64_t final_particle_count = num_particles_ + dim1_final_capacity + dim0_final_capacity;
        if (final_particle_count > total_capacity_)
        {
            resizeGlobalBuffers(final_particle_count + extra_capacity_);
        }

        auto particle_buffer_bytes = [](uint64_t capacity) -> size_t
        {
            return static_cast<size_t>(capacity) * (6 * sizeof(Real) + sizeof(long));
        };

        size_t workspace_bytes =
            particle_buffer_bytes(route1_capacity) +
            particle_buffer_bytes(route7_capacity) +
            particle_buffer_bytes(dim1_final_capacity);
        void * comm_workspace = nullptr;
        bool comm_workspace_private = false;
#ifndef DEBUG_RADIX_SORT
        cudaStream_t & comm_stream = pcl_stream;
#else
        cudaStream_t comm_stream = 0;
#endif

        nvtxRangePushA("moveParticles: CUDA-aware particle device workspace");
        acquireTemporaryWorkspace(&comm_workspace, &comm_workspace_private, workspace_bytes, comm_stream, "perfParticles::moveParticles");
        auto workspace_success = cudaStreamSynchronize(comm_stream);
        if (workspace_success != cudaSuccess)
        {
            nvtxRangePop();
            std::cerr << "CUDA workspace acquisition failed: " << cudaGetErrorString(workspace_success) << std::endl;
            throw std::runtime_error("Error acquiring CUDA-aware particle workspace");
        }
        nvtxRangePop();

        char * workspace_cursor = static_cast<char *>(comm_workspace);
        auto carve_particle_buffer = [&](device_particle_buffer & buffer, uint64_t capacity)
        {
            buffer.count = 0;
            buffer.capacity = capacity;
            if (capacity == 0)
            {
                buffer.p = nullptr;
                buffer.q = nullptr;
                buffer.other = nullptr;
                return;
            }

            buffer.p = reinterpret_cast<Real *>(workspace_cursor);
            workspace_cursor += 3 * capacity * sizeof(Real);
            buffer.q = reinterpret_cast<Real *>(workspace_cursor);
            workspace_cursor += 3 * capacity * sizeof(Real);
            buffer.other = reinterpret_cast<long *>(workspace_cursor);
            workspace_cursor += capacity * sizeof(long);
        };

        device_particle_buffer route1;
        device_particle_buffer route7;
        // First-hop sends still read the outgoing global slices during dim1 exchange.
        // Stage final dim1 arrivals until those sends have completed.
        device_particle_buffer dim1_final;
        carve_particle_buffer(route1, route1_capacity);
        carve_particle_buffer(route7, route7_capacity);
        carve_particle_buffer(dim1_final, dim1_final_capacity);

        auto copy_class_to_route = [&](device_particle_buffer & route, int proc)
        {
            uint64_t count = static_cast<uint64_t>(class_count(proc));
            if (count == 0) return;

            cudaMemcpy(route.p, p + 3 * send_begin[proc], 3 * count * sizeof(Real), cudaMemcpyDeviceToDevice);
            cudaMemcpy(route.q, q + 3 * send_begin[proc], 3 * count * sizeof(Real), cudaMemcpyDeviceToDevice);
            cudaMemcpy(route.other, other + send_begin[proc], count * sizeof(long), cudaMemcpyDeviceToDevice);
            route.count = count;
        };

        nvtxRangePushA("moveParticles: CUDA-aware particle route prep");
        copy_class_to_route(route1, 1);
        copy_class_to_route(route7, 7);
        for (int proc = 0; proc < 9; proc++)
        {
            extra_buffer_[proc].count = 0;
        }
        nvtxRangePop();

        auto send_dim1_class = [&](int proc, int to)
        {
            long count = class_count(proc);
            if (count == 0) return;

            parallel.send_dim1(p + 3 * send_begin[proc], 3 * count, to);
            parallel.send_dim1(q + 3 * send_begin[proc], 3 * count, to);
            parallel.send_dim1(other + send_begin[proc], count, to);
        };

        auto receive_dim1_buffer = [&](device_particle_buffer & buffer, long count, int from)
        {
            if (count == 0) return;

            parallel.receive_dim1(buffer.p + 3 * buffer.count, 3 * count, from);
            parallel.receive_dim1(buffer.q + 3 * buffer.count, 3 * count, from);
            parallel.receive_dim1(buffer.other + buffer.count, count, from);
            buffer.count += static_cast<uint64_t>(count);
        };

        auto receive_dim1_phase = [&](int phase, int from)
        {
            receive_dim1_buffer(route1, dim1_recv_sizes[phase][0], from);
            receive_dim1_buffer(dim1_final, dim1_recv_sizes[phase][1], from);
            receive_dim1_buffer(route7, dim1_recv_sizes[phase][2], from);
        };

        nvtxRangePushA("moveParticles: CUDA-aware particle payload MPI");

        if (parallel.grid_rank()[1] % 2 == 0)
        {
            int to = (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1];
            send_dim1_class(2, to);
            send_dim1_class(5, to);
            send_dim1_class(8, to);
        }
        if (parallel.grid_rank()[1] % 2 == 1 || (parallel.grid_rank()[1] == 0 && parallel.grid_size()[1] % 2 == 1))
        {
            receive_dim1_phase(0, (parallel.grid_rank()[1] + parallel.grid_size()[1] - 1) % parallel.grid_size()[1]);
        }

        if (parallel.grid_rank()[1] % 2 == 1)
        {
            int to = parallel.grid_rank()[1] - 1;
            send_dim1_class(0, to);
            send_dim1_class(3, to);
            send_dim1_class(6, to);
        }
        else
        {
            receive_dim1_phase(1, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
        }

        if (parallel.grid_rank()[1] % 2 == 0)
        {
            int to = (parallel.grid_rank()[1] + parallel.grid_size()[1] - 1) % parallel.grid_size()[1];
            send_dim1_class(0, to);
            send_dim1_class(3, to);
            send_dim1_class(6, to);
        }
        else
        {
            receive_dim1_phase(2, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
        }

        if (parallel.grid_rank()[1] % 2 == 1)
        {
            int to = (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1];
            send_dim1_class(2, to);
            send_dim1_class(5, to);
            send_dim1_class(8, to);
        }
        else if (parallel.grid_rank()[1] > 0 || parallel.grid_size()[1] % 2 == 0)
        {
            receive_dim1_phase(3, (parallel.grid_rank()[1] + parallel.grid_size()[1] - 1) % parallel.grid_size()[1]);
        }

        if (dim1_final.count > 0)
        {
            nvtxRangePushA("moveParticles: CUDA-aware particle append dim1 arrivals");
            cudaMemcpy(p + 3 * num_particles_, dim1_final.p, 3 * dim1_final.count * sizeof(Real), cudaMemcpyDeviceToDevice);
            cudaMemcpy(q + 3 * num_particles_, dim1_final.q, 3 * dim1_final.count * sizeof(Real), cudaMemcpyDeviceToDevice);
            cudaMemcpy(other + num_particles_, dim1_final.other, dim1_final.count * sizeof(long), cudaMemcpyDeviceToDevice);
            num_particles_ += dim1_final.count;
            nvtxRangePop();
        }

        auto send_dim0_route = [&](device_particle_buffer & route, int to)
        {
            if (route.count == 0) return;

            parallel.send_dim0(route.p, 3 * route.count, to);
            parallel.send_dim0(route.q, 3 * route.count, to);
            parallel.send_dim0(route.other, route.count, to);
        };

        auto receive_dim0_final = [&](long count, int from)
        {
            if (count == 0) return;

            parallel.receive_dim0(p + 3 * num_particles_, 3 * count, from);
            parallel.receive_dim0(q + 3 * num_particles_, 3 * count, from);
            parallel.receive_dim0(other + num_particles_, count, from);
            num_particles_ += static_cast<uint64_t>(count);
        };

        if (parallel.grid_rank()[0] % 2 == 0)
        {
            send_dim0_route(route7, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
        }
        if (parallel.grid_rank()[0] % 2 == 1 || (parallel.grid_rank()[0] == 0 && parallel.grid_size()[0] % 2 == 1))
        {
            receive_dim0_final(dim0_recv_sizes[0], (parallel.grid_rank()[0] + parallel.grid_size()[0] - 1) % parallel.grid_size()[0]);
        }

        if (parallel.grid_rank()[0] % 2 == 1)
        {
            send_dim0_route(route1, parallel.grid_rank()[0] - 1);
        }
        else
        {
            receive_dim0_final(dim0_recv_sizes[1], (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
        }

        if (parallel.grid_rank()[0] % 2 == 0)
        {
            send_dim0_route(route1, (parallel.grid_rank()[0] + parallel.grid_size()[0] - 1) % parallel.grid_size()[0]);
        }
        else
        {
            receive_dim0_final(dim0_recv_sizes[2], (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
        }

        if (parallel.grid_rank()[0] % 2 == 1)
        {
            send_dim0_route(route7, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
        }
        else if (parallel.grid_rank()[0] > 0 || parallel.grid_size()[0] % 2 == 0)
        {
            receive_dim0_final(dim0_recv_sizes[3], (parallel.grid_rank()[0] + parallel.grid_size()[0] - 1) % parallel.grid_size()[0]);
        }

        nvtxRangePop();

        releaseTemporaryWorkspace(&comm_workspace, &comm_workspace_private, comm_stream);
    }
    else
    {
        nvtxRangePushA("moveParticles: host-staged particle migration");
        nvtxRangePushA("moveParticles: set up communication buffers");

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
    nvtxRangePop();
    nvtxRangePop();
    }

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

    computeSortIndices(d_keys, d_indices_in, &d_temp, &d_temp_private, pcl_stream);

    nvtxRangePop();

    reorderParticles(d_indices_in, d_temp, pcl_stream);

    releaseTemporaryWorkspace(&d_temp, &d_temp_private, pcl_stream);
    cudaFreeAsync(d_indices_in, pcl_stream);

    update_pointers<<<(num_row_buffers_+137)/128, 128, 0, pcl_stream>>>(this, nullptr);

    success = cudaStreamSynchronize(pcl_stream);
    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error after cudaStreamSynchronize in moveParticles");
    }

    cudaStreamDestroy(pcl_stream);
#endif

    nvtxRangePop(); // moveParticles
}

// project particles
template <typename part, typename part_info>
__global__ void project_particles(perfParticles<part, part_info> * pcl, Real * target, int projection_order, long start_idx, long jump0, long jump1, long jump2, int stencil_k)
{
    int row = blockIdx.x;
    int thread_id = threadIdx.x;

    for (int idx = thread_id; idx < pcl->row_buffers_[row].count; idx += 128)
    {
        pcl->project_particle(target, projection_order, start_idx, jump0, jump1, jump2, row, stencil_k, idx);
    }
}

// density projection
template <typename part, typename part_info>
void perfParticles<part, part_info>::meshprojection_project(Field<Real> * target, int projection_order)
{
    long start_idx = target->lattice().siteFirst();
    long jump0, jump1, jump2;

    Real * data = target->data();

    jump0 = target->lattice().jump(0);
    jump1 = target->lattice().jump(1);
    jump2 = target->lattice().jump(2);

    // Loop over the projection stencil
    for (int k = 0; k < 4; k++)
    {
        project_particles<<<num_row_buffers_, 128>>>(this, data, projection_order, start_idx, jump0, jump1, jump2, k);
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
__device__ void perfParticles<part, part_info>::project_particle(Real * target, int projection_order, long start_idx, long jump0, long jump1, long jump2, int row, int stencil_k, int idx)
{
    int coord[3];
    Real frac[3];
    Real dx = boxSize_[0] / lat_size_[0];

    Real weight = part_global_info_.mass / (dx * dx * dx);

    coord[1] = static_cast<int>(floor(row_buffers_[row].p[1]*lat_size_[1]));
    coord[2] = static_cast<int>(floor(row_buffers_[row].p[2]*lat_size_[2]));
    coord[0] = 0;

    long site = start_idx + coord[0]*jump0 + (coord[1] - coordSkip_[1])*jump1 + (coord[2] - coordSkip_[0])*jump2;

    if (stencil_k % 2)
    {
        site += jump1;
        coord[1]++;
    }
    if (stencil_k / 2)
    {
        site += jump2;
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
