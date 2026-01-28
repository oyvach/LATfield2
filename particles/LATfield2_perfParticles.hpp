#ifndef LATFIELD2_PERFPARTICLES_HPP
#define LATFIELD2_PERFPARTICLES_HPP

#include <algorithm>

#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>

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
    // Reallocate the arrays using cuda_reallocManaged
 /*   p_type * new_p = cuda_reallocManaged<p_type>(p, capacity * 3, new_capacity * 3);
    q_type * new_q = cuda_reallocManaged<q_type>(q, capacity * 3, new_capacity * 3);
    other_type * new_other = cuda_reallocManaged<other_type>(other, capacity, new_capacity);
    auto success = cudaDeviceSynchronize();

    if (success != cudaSuccess)
    {
        std::cerr << "CUDA memory reallocation failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in cuda_reallocManaged");
    }

    // Update the pointers
    p = new_p;
    q = new_q;
    other = new_other;
    // Update capacity
    capacity = new_capacity;*/
    resize(new_capacity);
}

template <typename part, typename part_info>
class perfParticles;

template <typename part, typename part_info>
__global__ void update_pointers(perfParticles<part, part_info> * pcl, unsigned long long int * send_begin);

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
        //row_buffer<Real, Real, long> add_buffer_;
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

            /*total_capacity_ = new_total_capacity;

            Real * new_p = (Real *) realloc(p, total_capacity_ * 3 * sizeof(Real));
            Real * new_q = (Real *) realloc(q, total_capacity_ * 3 * sizeof(Real));
            long * new_other = (long *) realloc(other, total_capacity_ * sizeof(long));

            if (new_p == NULL || new_q == NULL || new_other == NULL)
            {
                throw std::runtime_error("Memory allocation failed");
            }

            p = new_p;
            q = new_q;
            other = new_other;*/

            p = cuda_realloc(p, total_capacity_ * 3, new_total_capacity * 3);
            q = cuda_realloc(q, total_capacity_ * 3, new_total_capacity * 3);
            other = cuda_realloc(other, total_capacity_, new_total_capacity);

            total_capacity_ = new_total_capacity;
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

        template <typename part2, typename part_info2, typename UpdateFunct>
        friend __global__ void update_particles(perfParticles<part2, part_info2> * pcl, UpdateFunct update_funct, double dtau, Field<Real> ** fields, int nfields, double * params, double * output, int * reduce_type, int noutput, Real * v2, bool copyout);

        template <typename part2, typename part_info2>
        friend __global__ void project_particles(perfParticles<part2, part_info2> * pcl, Real * target, int projection_order, long start_idx, const long * jump, int stencil_k);
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
    /*add_buffer_.p = nullptr;
    add_buffer_.q = nullptr;
    add_buffer_.other = nullptr;
    add_buffer_.count = 0;
    add_buffer_.capacity = 0;*/
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
        /*free(row_buffers_);
        free(p);
        free(q);
        free(other);*/
        cudaFree(row_buffers_);
        cudaFree(p);
        cudaFree(q);
        cudaFree(other);
    }

    /*cudaFree(add_buffer_.p);
    cudaFree(add_buffer_.q);
    cudaFree(add_buffer_.other);*/
    for (int i = 0; i < 9; i++)
    {
        //cudaFree(extra_buffer_[i].p);
        //cudaFree(extra_buffer_[i].q);
        //cudaFree(extra_buffer_[i].other);
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
    //row_buffers_ = (row_buffer<Real, Real, long> *)malloc(num_row_buffers_ * sizeof(row_buffer<Real, Real, long>));
    auto success = cudaMalloc(&row_buffers_, num_row_buffers_ * sizeof(row_buffer<Real, Real, long>));
    // Check if the allocation was successful
    if (success != cudaSuccess || row_buffers_ == nullptr)
    {
        std::cerr << " proc#" << parallel.rank() << " cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
        // If not, throw an exception
        throw std::runtime_error("Memory allocation failed for row_buffers_ in perfParticles::initialize");
    }
    // Initialize the row buffers
    //for (uint32_t i = 0; i < num_row_buffers_; i++)
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
    //p = (Real *)malloc(total_capacity_ * 3 * sizeof(Real));
    //q = (Real *)malloc(total_capacity_ * 3 * sizeof(Real));
    //other = (long *)malloc(total_capacity_ * sizeof(long));
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
        //extra_buffer_[i].capacity = ((i == 4 || i % 2) ? extra_capacity_ : (extra_capacity_)); // corners need less capacity
        extra_buffer_[i].sorted = true;

        extra_buffer_[i].p = (Real *) malloc(extra_buffer_[i].capacity * 3 * sizeof(Real));
        extra_buffer_[i].q = (Real *) malloc(extra_buffer_[i].capacity * 3 * sizeof(Real));
        extra_buffer_[i].other = (long *) malloc(extra_buffer_[i].capacity * sizeof(long));

        /*success = cudaMallocManaged(&extra_buffer_[i].p, extra_buffer_[i].capacity * 3 * sizeof(Real));

        if (success != cudaSuccess)
        {
            std::cerr << "cudaMallocManaged failed: " << cudaGetErrorString(success) << std::endl;
            throw std::runtime_error("Memory allocation failed in perfParticles::initialize");
        }

        success = cudaMallocManaged(&extra_buffer_[i].q, extra_buffer_[i].capacity * 3 * sizeof(Real));

        if (success != cudaSuccess)
        {
            std::cerr << "cudaMallocManaged failed: " << cudaGetErrorString(success) << std::endl;
            throw std::runtime_error("Memory allocation failed in perfParticles::initialize");
        }

        success = cudaMallocManaged(&extra_buffer_[i].other, extra_buffer_[i].capacity * sizeof(long));

        if (success != cudaSuccess)
        {
            std::cerr << "cudaMallocManaged failed: " << cudaGetErrorString(success) << std::endl;
            throw std::runtime_error("Memory allocation failed in perfParticles::initialize");
        }*/
    }
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

    // check if we have space to add particle
    /*if (num_particles_ >= total_capacity_)
    {
        // realloc
        //total_capacity_ += extra_capacity_;
        //p = (Real *)realloc(p, total_capacity_ * 3 * sizeof(Real));
        //q = (Real *)realloc(q, total_capacity_ * 3 * sizeof(Real));
        //other = (long *)realloc(other, total_capacity_ * sizeof(long));
        resizeGlobalBuffers(total_capacity_ + extra_capacity_);
    }*/

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

    // Add the particle to the global array
    /*for (int i = 0; i < 3; i++)
    {
        p[3 * num_particles_ + i] = newPart.pos[i];
        q[3 * num_particles_ + i] = newPart.vel[i];
    }
    other[num_particles_] = newPart.ID;*/

    //cudaMemcpy(p + 3 * num_particles_, newPart.pos, 3 * sizeof(Real), cudaMemcpyDefault);
    //cudaMemcpy(q + 3 * num_particles_, newPart.vel, 3 * sizeof(Real), cudaMemcpyDefault);
    //cudaMemcpy(other + num_particles_, &newPart.ID, sizeof(long), cudaMemcpyDefault);

    // Increment the particle count
    // num_particles_++;

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

// update row buffers
template <typename part, typename part_info>
void perfParticles<part, part_info>::updateRowBuffers(unsigned long long int * send_begin)
{
    // check if adding buffer is empty, otherwise add particles to global buffer
    nvtxRangePushA("updateRowBuffers: add particles");

    flushExtraBuffer(4);

    nvtxRangePop();

    // tuple-pointer to particle properties
    tripleReal * pos = (tripleReal *) p;
    tripleReal * vel = (tripleReal *) q;

    // create zip iterator from tuples
    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(pos, vel, other));
    auto zip_end = thrust::make_zip_iterator(thrust::make_tuple(pos + num_particles_, vel + num_particles_, other + num_particles_));

    using zip_iter = typename decltype(zip_begin)::value_type;

    nvtxRangePushA("updateRowBuffers: sort");

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

    nvtxRangePop();

    nvtxRangePushA("updateRowBuffers: update pointers");

    // loop over row buffers
    /*for (int row = 0; row < num_row_buffers_; row++)
    {
        auto row_begin = thrust::lower_bound(thrust::device, zip_begin, zip_end, row, [this] __device__ (zip_iter a, int row) {
            auto pos_a = thrust::get<0>(a);
            return computeRow(((Real *) (&pos_a))[1], ((Real *) (&pos_a))[2]) < row;
        });

        auto row_end = thrust::upper_bound(thrust::device, row_begin, zip_end, row, [this] __device__ (int row, zip_iter a) {
            auto pos_a = thrust::get<0>(a);
            return row < computeRow(((Real *) (&pos_a))[1], ((Real *) (&pos_a))[2]);
        });

        row_buffers_[row].count = row_end - row_begin;

        if (row_buffers_[row].count > 0)
        {
            row_buffers_[row].p = (Real *) & thrust::get<0>(*row_begin);
            row_buffers_[row].q = (Real *) & thrust::get<1>(*row_begin);
            row_buffers_[row].other = & thrust::get<2>(*row_begin);
        }
        else
        {
            row_buffers_[row].p = nullptr;
            row_buffers_[row].q = nullptr;
            row_buffers_[row].other = nullptr;
        }

        row_buffers_[row].sorted = true;
    }*/

    if (send_begin != nullptr)
    {
        update_pointers<<<(num_row_buffers_+137)/128, 128>>>(this, send_begin);
    }
    else
    {
        update_pointers<<<(num_row_buffers_+127)/128, 128>>>(this, nullptr);
    }

    auto success = cudaDeviceSynchronize();

    if (success != cudaSuccess)
    {
        std::cerr << "CUDA kernel failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Error in CUDA kernel update_pointers");
    }

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

    update_particles<<<num_row_buffers_, 128>>>(this, update_funct, dtau, fields, nfields, params, output, reduce_type, noutput, &v2, copyout);

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
    nvtxRangePushA("moveParticles: updateParticles");

    updateParticles(move_funct, dtau, fields, nfields, params, output, reduce_type, noutput);

    nvtxRangePop();

    /*for (int i = 0; i < num_row_buffers_; i++)
    {
        row_buffers_[i].sorted = false;
    }*/

    thrust::for_each(thrust::device, row_buffers_, row_buffers_ + num_row_buffers_, [] __device__ (row_buffer<Real, Real, long> & rb) { rb.sorted = false;});

    rows_sorted_ = false;

    nvtxRangePushA("moveParticles: updateRowBuffers");

    unsigned long long int send_begin[10];

    updateRowBuffers(&send_begin[0]);

    nvtxRangePop();

    // memory movement

    //row_buffer<Real, Real, long> send_buffer[9];

    nvtxRangePushA("moveParticles: set up communication buffers");

    /*tripleReal * pos = (tripleReal *) p;
    tripleReal * vel = (tripleReal *) q;

    // create zip iterator from tuples
    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(pos, vel, other));
    auto zip_end = thrust::make_zip_iterator(thrust::make_tuple(pos + num_particles_, vel + num_particles_, other + num_particles_));

    using zip_iter = typename decltype(zip_begin)::value_type;

    auto rec_begin = thrust::lower_bound(thrust::device, zip_begin, zip_end, num_row_buffers_, [this] __host__ __device__ (zip_iter a, int row) {
        auto pos_a = thrust::get<0>(a);
        return computeRow(((Real *) (&pos_a))[1], ((Real *) (&pos_a))[2]) < row;
    });*/

    /*send_buffer[4].count = 0;
    send_buffer[4].capacity = -1;

    send_buffer[4].p = (Real *) & thrust::get<0>(*rec_begin);
    send_buffer[4].q = (Real *) & thrust::get<1>(*rec_begin);
    send_buffer[4].other = & thrust::get<2>(*rec_begin);*/

    num_particles_ = send_begin[0]; //rec_begin - zip_begin;

    // loop over send buffers
    for (int proc = 0; proc < 9; proc++)
    {
        //if (proc == 4)
        //{
        //    continue;
        //}
        //else
        //{
        /*auto send_begin = thrust::lower_bound(thrust::device, rec_begin, zip_end, num_row_buffers_+proc, [this] __host__ __device__ (zip_iter a, int row) {
            auto pos_a = thrust::get<0>(a);
            return computeRow(((Real *) (&pos_a))[1], ((Real *) (&pos_a))[2]) < row;
        });

        auto send_end = thrust::upper_bound(thrust::device, send_begin, zip_end, num_row_buffers_+proc, [this] __host__ __device__ (int row, zip_iter a) {
            auto pos_a = thrust::get<0>(a);
            return row < computeRow(((Real *) (&pos_a))[1], ((Real *) (&pos_a))[2]);
        });*/

        extra_buffer_[proc].count = send_begin[proc+1] - send_begin[proc]; //send_end - send_begin;

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

        //extra_buffer_[proc].capacity = extra_buffer_[proc].count;

        /*if (extra_buffer_[proc].count > 0)
        {
            // use cudaMallocManaged
            cudaMallocManaged(&send_buffer[proc].p, send_buffer[proc].count * 3 * sizeof(Real));

            send_buffer[proc].p = (Real *)malloc(send_buffer[proc].count * 3 * sizeof(Real));
            send_buffer[proc].q = (Real *)malloc(send_buffer[proc].count * 3 * sizeof(Real));
            send_buffer[proc].other = (long *)malloc(send_buffer[proc].count * sizeof(long));


            // copy particles to send buffer using memcpy
            //memcpy(send_buffer[proc].p, & thrust::get<0>(*send_begin), send_buffer[proc].count * 3 * sizeof(Real));
            //memcpy(send_buffer[proc].q, & thrust::get<1>(*send_begin), send_buffer[proc].count * 3 * sizeof(Real));
            //memcpy(send_buffer[proc].other, & thrust::get<2>(*send_begin), send_buffer[proc].count * sizeof(long));

            cudaMemcpy(send_buffer[proc].p, & thrust::get<0>(*send_begin), send_buffer[proc].count * 3 * sizeof(Real), cudaMemcpyDefault);
            cudaMemcpy(send_buffer[proc].q, & thrust::get<1>(*send_begin), send_buffer[proc].count * 3 * sizeof(Real), cudaMemcpyDefault);
            cudaMemcpy(send_buffer[proc].other, & thrust::get<2>(*send_begin), send_buffer[proc].count * sizeof(long), cudaMemcpyDefault);
        }
        else
        {
            send_buffer[proc].p = nullptr;
            send_buffer[proc].q = nullptr;
            send_buffer[proc].other = nullptr;
        }*/
        //}
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

        // free send buffer particle arrays
        /*free(send_buffer[2].p);
        free(send_buffer[5].p);
        free(send_buffer[8].p);
        free(send_buffer[2].q);
        free(send_buffer[5].q);
        free(send_buffer[8].q);
        free(send_buffer[2].other);
        free(send_buffer[5].other);
        free(send_buffer[8].other);*/
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
            //extra_buffer_[1].resize(buffer_sizes[0]+extra_buffer_[1].count);
            parallel.receive_dim1(extra_buffer_[1].p+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].q+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].other+extra_buffer_[1].count, buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            extra_buffer_[1].count += buffer_sizes[0];
        }
        if (buffer_sizes[1] > 0)
        {
            /*if (num_particles_ + buffer_sizes[1] > total_capacity_)
            {
                resizeGlobalBuffers(total_capacity_ + buffer_sizes[1] + extra_capacity_);
                send_buffer[4].p = p + 3*num_particles_;
                send_buffer[4].q = q + 3*num_particles_;
                send_buffer[4].other = other + num_particles_;
            }*/
           
            if (extra_buffer_[4].count + buffer_sizes[1] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[1] + extra_capacity_);
            }
            //extra_buffer_[4].resize(buffer_sizes[1]+extra_buffer_[4].count);
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
            //extra_buffer_[7].resize(buffer_sizes[2]+extra_buffer_[7].count);
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

        // free send buffer particle arrays
        /*free(send_buffer[0].p);
        free(send_buffer[3].p);
        free(send_buffer[6].p);
        free(send_buffer[0].q);
        free(send_buffer[3].q);
        free(send_buffer[6].q);
        free(send_buffer[0].other);
        free(send_buffer[3].other);
        free(send_buffer[6].other);*/
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
            //extra_buffer_[1].resize(buffer_sizes[0]+extra_buffer_[1].count);
            parallel.receive_dim1(extra_buffer_[1].p+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].q+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].other+extra_buffer_[1].count, buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            extra_buffer_[1].count += buffer_sizes[0];
        }
        if (buffer_sizes[1] > 0)
        {
            /*if (num_particles_ + send_buffer[4].count + buffer_sizes[1] > total_capacity_)
            {
                resizeGlobalBuffers(total_capacity_ + send_buffer[4].count + buffer_sizes[1] + extra_capacity_);
                send_buffer[4].p = p + 3*num_particles_;
                send_buffer[4].q = q + 3*num_particles_;
                send_buffer[4].other = other + num_particles_;
            }*/

            if (extra_buffer_[4].count + buffer_sizes[1] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[1] + extra_capacity_);
            }
            //extra_buffer_[4].resize(buffer_sizes[1]+extra_buffer_[4].count);
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
            //extra_buffer_[7].resize(buffer_sizes[2]+extra_buffer_[7].count);
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

        // free send buffer particle arrays
        /*free(send_buffer[0].p);
        free(send_buffer[3].p);
        free(send_buffer[6].p);
        free(send_buffer[0].q);
        free(send_buffer[3].q);
        free(send_buffer[6].q);
        free(send_buffer[0].other);
        free(send_buffer[3].other);
        free(send_buffer[6].other);*/
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
            //extra_buffer_[1].resize(buffer_sizes[0]+extra_buffer_[1].count);
            parallel.receive_dim1(extra_buffer_[1].p+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].q+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].other+extra_buffer_[1].count, buffer_sizes[0], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
            extra_buffer_[1].count += buffer_sizes[0];
        }
        if (buffer_sizes[1] > 0)
        {
            /*if (num_particles_ + send_buffer[4].count + buffer_sizes[1] > total_capacity_)
            {
                resizeGlobalBuffers(total_capacity_ + send_buffer[4].count + buffer_sizes[1] + extra_capacity_);
                send_buffer[4].p = p + 3*num_particles_;
                send_buffer[4].q = q + 3*num_particles_;
                send_buffer[4].other = other + num_particles_;
            }*/

            if (extra_buffer_[4].count + buffer_sizes[1] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[1] + extra_capacity_);
            }
            //extra_buffer_[4].resize(buffer_sizes[1]+extra_buffer_[4].count);
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
            //extra_buffer_[7].resize(buffer_sizes[2]+extra_buffer_[7].count);
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

        // free send buffer particle arrays
        /*free(send_buffer[2].p);
        free(send_buffer[5].p);
        free(send_buffer[8].p);
        free(send_buffer[2].q);
        free(send_buffer[5].q);
        free(send_buffer[8].q);
        free(send_buffer[2].other);
        free(send_buffer[5].other);
        free(send_buffer[8].other);*/
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
            //extra_buffer_[1].resize(buffer_sizes[0]+extra_buffer_[1].count);
            parallel.receive_dim1(extra_buffer_[1].p+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].q+3*extra_buffer_[1].count, 3*buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            parallel.receive_dim1(extra_buffer_[1].other+extra_buffer_[1].count, buffer_sizes[0], (parallel.grid_rank()[1]+parallel.grid_size()[1]-1) % parallel.grid_size()[1]);
            extra_buffer_[1].count += buffer_sizes[0];
        }
        if (buffer_sizes[1] > 0)
        {
            /*if (num_particles_ + send_buffer[4].count + buffer_sizes[1] > total_capacity_)
            {
                resizeGlobalBuffers(total_capacity_ + send_buffer[4].count + buffer_sizes[1] + extra_capacity_);
                send_buffer[4].p = p + 3*num_particles_;
                send_buffer[4].q = q + 3*num_particles_;
                send_buffer[4].other = other + num_particles_;
            }*/

            if (extra_buffer_[4].count + buffer_sizes[1] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[1] + extra_capacity_);
            }
            //extra_buffer_[4].resize(buffer_sizes[1]+extra_buffer_[4].count);
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
            //extra_buffer_[7].resize(buffer_sizes[2]+extra_buffer_[7].count);
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
            //num_particles_ -= buffer_sizes[0];
            parallel.send_dim0(extra_buffer_[7].p, 3*extra_buffer_[7].count, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.send_dim0(extra_buffer_[7].q, 3*extra_buffer_[7].count, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
            parallel.send_dim0(extra_buffer_[7].other, extra_buffer_[7].count, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
        }

        // free send buffer particle arrays
        /*free(send_buffer[7].p);
        free(send_buffer[7].q);
        free(send_buffer[7].other);*/
    }
    
    if (parallel.grid_rank()[0] % 2 == 1 || (parallel.grid_rank()[0] == 0 && parallel.grid_size()[0] % 2 == 1)) // odd rank or rank 0 with odd grid size
    {
        parallel.receive_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            /*if (num_particles_ + send_buffer[4].count + buffer_sizes[0] > total_capacity_)
            {
                resizeGlobalBuffers(total_capacity_ + send_buffer[4].count + buffer_sizes[0] + extra_capacity_);
                send_buffer[4].p = p + 3*num_particles_;
                send_buffer[4].q = q + 3*num_particles_;
                send_buffer[4].other = other + num_particles_;
            }*/

            if (extra_buffer_[4].count + buffer_sizes[0] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[0] + extra_capacity_);
            }
            //extra_buffer_[4].resize(buffer_sizes[0]+extra_buffer_[4].count);
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

        // free send buffer particle arrays
        /*free(send_buffer[1].p);
        free(send_buffer[1].q);
        free(send_buffer[1].other);*/
    }
    else
    {
        parallel.receive_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            /*if (num_particles_ + send_buffer[4].count + buffer_sizes[0] > total_capacity_)
            {
                resizeGlobalBuffers(total_capacity_ + send_buffer[4].count + buffer_sizes[0] + extra_capacity_);
                send_buffer[4].p = p + 3*num_particles_;
                send_buffer[4].q = q + 3*num_particles_;
                send_buffer[4].other = other + num_particles_;
            }*/

            if (extra_buffer_[4].count + buffer_sizes[0] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[0] + extra_capacity_);
            }
            //extra_buffer_[4].resize(buffer_sizes[0]+extra_buffer_[4].count);
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

        // free send buffer particle arrays
        /*free(send_buffer[1].p);
        free(send_buffer[1].q);
        free(send_buffer[1].other);*/
    }
    else
    {
        parallel.receive_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            /*if (num_particles_ + send_buffer[4].count + buffer_sizes[0] > total_capacity_)
            {
                resizeGlobalBuffers(total_capacity_ + send_buffer[4].count + buffer_sizes[0] + extra_capacity_);
                send_buffer[4].p = p + 3*num_particles_;
                send_buffer[4].q = q + 3*num_particles_;
                send_buffer[4].other = other + num_particles_;
            }*/

            if (extra_buffer_[4].count + buffer_sizes[0] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[0] + extra_capacity_);
            }
            //extra_buffer_[4].resize(buffer_sizes[0]+extra_buffer_[4].count);
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

        // free send buffer particle arrays
        /*free(send_buffer[7].p);
        free(send_buffer[7].q);
        free(send_buffer[7].other);*/
    }
    else if (parallel.grid_rank()[0] > 0 || parallel.grid_size()[0] % 2 == 0)
    {
        parallel.receive_dim0(buffer_sizes, 1, (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);

        if (buffer_sizes[0] > 0)
        {
            /*if (num_particles_ + send_buffer[4].count + buffer_sizes[0] > total_capacity_)
            {
                resizeGlobalBuffers(total_capacity_ + send_buffer[4].count + buffer_sizes[0] + extra_capacity_);
                send_buffer[4].p = p + 3*num_particles_;
                send_buffer[4].q = q + 3*num_particles_;
                send_buffer[4].other = other + num_particles_;
            }*/

            if (extra_buffer_[4].count + buffer_sizes[0] > extra_buffer_[4].capacity)
            {
                extra_buffer_[4].resizeManaged(extra_buffer_[4].count + buffer_sizes[0] + extra_capacity_);
            }
            //extra_buffer_[4].resize(buffer_sizes[0]+extra_buffer_[4].count);
            parallel.receive_dim0(extra_buffer_[4].p+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].q+3*extra_buffer_[4].count, 3*buffer_sizes[0], (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            parallel.receive_dim0(extra_buffer_[4].other+extra_buffer_[4].count, buffer_sizes[0], (parallel.grid_rank()[0]+parallel.grid_size()[0]-1) % parallel.grid_size()[0]);
            extra_buffer_[4].count += buffer_sizes[0];
        }
    }

    nvtxRangePop();

    // ingest particles from send_buffer[4]

    nvtxRangePushA("moveParticles: ingest received particles");

    /*// check if there is enough capacity in the global buffers
    if (num_particles_ + send_buffer[4].count > total_capacity_)
    {
        nvtxRangePushA("moveParticles: resize global buffers");
        resizeGlobalBuffers(total_capacity_ + send_buffer[4].count + extra_capacity_);
        nvtxRangePop();
    }

    nvtxRangePushA("moveParticles: copy particles to global buffers");

    // copy particles from send_buffer[4] to global buffers
    cudaMemcpy(p+3*num_particles_, send_buffer[4].p, 3*send_buffer[4].count*sizeof(Real), cudaMemcpyDefault);
    cudaMemcpy(q+3*num_particles_, send_buffer[4].q, 3*send_buffer[4].count*sizeof(Real), cudaMemcpyDefault);
    cudaMemcpy(other+num_particles_, send_buffer[4].other, send_buffer[4].count*sizeof(long), cudaMemcpyDefault);

    nvtxRangePop();

    // free send buffer particle arrays
    free(send_buffer[4].p);
    free(send_buffer[4].q);
    free(send_buffer[4].other);

    num_particles_ += send_buffer[4].count;*/

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

    nvtxRangePushA("moveParticles: update row buffers");

    updateRowBuffers();

    nvtxRangePop();

    /*// check consistency of particle counts
    uint64_t count = 0;

    for (int row = 0; row < num_row_buffers_; row++)
    {
        count += row_buffers_[row].count;
    }

    if (count != num_particles_)
    {
        cerr << " proc#" << parallel.rank() << ": Inconsistent particle counts after MPI communication: " << count << " (sum over row buffers) != " << num_particles_  << " (num_particles_); position of first particle not assigned a row: (" << p[3*count] << ", " << p[3*count+1] << ", " << p[3*count+2] << ") which is in fictitious row #" << computeRow(p[3*count+1], p[3*count+2])-num_row_buffers_ << endl;
        throw std::runtime_error("Inconsistent particle counts");
    }*/

    nvtxRangePop();
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
