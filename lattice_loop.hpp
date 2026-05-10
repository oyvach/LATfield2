#ifndef LATTICE_LOOP_HPP
#define LATTICE_LOOP_HPP

#include <cuda/atomic>

using namespace LATfield2;

// generic loop for lattice
template <typename ForEachFunct, int noutput = 0>
__global__ void lattice_for_each(ForEachFunct funct, int numpts, Field<Real> ** fields, int nfields, double * params,
                                double * output, int * reduce_type,void ** vparams = nullptr, int intoHalo = 0)
{
    int coord1 = blockIdx.x;
    int coord2 = blockIdx.y;
    int thread_id = threadIdx.x;

    alignas(8) double output_thread[(noutput > 0) ? noutput : 1];
    alignas(8) double output_site  [(noutput > 0) ? noutput : 1];

    Site * sites = nullptr;

    if (nfields > 0)
    {
        sites = (Site *)alloca(nfields * sizeof(Site));
    }

    if constexpr (noutput > 0)
    {
        #pragma unroll
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
    }

    for (int idx = thread_id; idx < numpts; idx += 128)
    {
        for (int i = 0; i < nfields; i++)
        {
            const long base = fields[i]->lattice().siteFirst() -
                              (long) intoHalo * (fields[i]->lattice().jump(0) + fields[i]->lattice().jump(1) + fields[i]->lattice().jump(2));

            sites[i] = Site(fields[i]->lattice(), base 
                            + idx*fields[i]->lattice().jump(0) 
                            + coord1*fields[i]->lattice().jump(1)
                            + coord2*fields[i]->lattice().jump(2));
        }
        funct(fields, sites, nfields, params, output_site, vparams);

        if constexpr (noutput > 0)
        {
            #pragma unroll
            for (int i = 0; i < noutput; i++)
            {
                if (reduce_type[i] & (SUM | SUM_LOCAL))
                {
                    output_thread[i] += output_site[i];
                }
                else if (reduce_type[i] & (MAX | MAX_LOCAL) && output_site[i] > output_thread[i])
                {
                    output_thread[i] = output_site[i];
                }
                else if (reduce_type[i] & (MIN | MIN_LOCAL) && output_site[i] < output_thread[i])
                {
                    output_thread[i] = output_site[i];
                }
            }
        }
    }

    if constexpr (noutput > 0)
    {
        __shared__ double smem[4 * noutput];

        constexpr int W = 32;
        const int lane   = threadIdx.x & (W-1);
        const int warpId = threadIdx.x / W;
        constexpr int nWarps = 4; // blockDim.x = 128

        unsigned mask = __activemask();

        double warp_val[noutput];

        #pragma unroll
        for (int i = 0; i < noutput; i++)
        {
            double v = output_thread[i];

            if (reduce_type[i] & SUM)
            {
                for (int ofs = 16; ofs > 0; ofs >>= 1)  // 16,8,4,2,1 = log2(32) steps
                    v += __shfl_down_sync(mask, v, ofs);
            }
            else if (reduce_type[i] & MAX)
            {
                for (int ofs = 16; ofs > 0; ofs >>= 1)  // 16,8,4,2,1 = log2(32) steps
                    v = fmax(v, __shfl_down_sync(mask, v, ofs));
            }
            else if (reduce_type[i] & MIN)
            {
                for (int ofs = 16; ofs > 0; ofs >>= 1)  // 16,8,4,2,1 = log2(32) steps
                    v = fmin(v, __shfl_down_sync(mask, v, ofs));
            }

            warp_val[i] = v;
        }

        if (lane == 0) {
            #pragma unroll
            for (int i = 0; i < noutput; ++i)
                smem[warpId * noutput + i] = warp_val[i];
        }
        __syncthreads();

        if (warpId == 0)
        {
            #pragma unroll
            for (int i = 0; i < noutput; i++)
            {
                double v = (lane < nWarps) ? smem[lane * noutput + i] : 0.0;

                if (reduce_type[i] & SUM)
                {
                    #pragma unroll
                    for (int ofs = 2; ofs > 0; ofs >>= 1)  // 2,1 = log2(4) steps, blockDim.x = 128 => nWarps = 4
                        v += __shfl_down_sync(mask, v, ofs);
                }
                else if (reduce_type[i] & MAX)
                {
                    #pragma unroll
                    for (int ofs = 2; ofs > 0; ofs >>= 1)  // 2,1 = log2(4) steps, blockDim.x = 128 => nWarps = 4
                        v = fmax(v, __shfl_down_sync(mask, v, ofs));
                }
                else if (reduce_type[i] & MIN)
                {
                    #pragma unroll
                    for (int ofs = 2; ofs > 0; ofs >>= 1)  // 2,1 = log2(4) steps, blockDim.x = 128 => nWarps = 4
                        v = fmin(v, __shfl_down_sync(mask, v, ofs));
                }

                if (lane == 0)
                {
                    cuda::atomic_ref<double, cuda::thread_scope_device> output_ref(output[i]);
                    if (reduce_type[i] & SUM)
                    {
                        output_ref.fetch_add(v);
                    }
                    else if (reduce_type[i] & MAX)
                    {
                        output_ref.fetch_max(v);
                    }
                    else if (reduce_type[i] & MIN)
                    {
                        output_ref.fetch_min(v);
                    }
                }
            }
        }
    }
}

template <typename ForEachHaloFunct>
__global__ void lattice_for_each_halo(ForEachHaloFunct funct, int numpts, int halo, int size0, int size1, int size2, Field<Real> ** fields, int nfields, double * params)
{
    int coord1 = blockIdx.x;
    int coord2 = blockIdx.y;
    int thread_id = threadIdx.x;

    Site * sites = nullptr;

    if (nfields > 0)
    {
        sites = (Site *)alloca(nfields * sizeof(Site));
    }

    for (int idx = thread_id; idx < numpts; idx += 128)
    {
        const bool in_halo = (idx < halo || idx >= size0 + halo ||
                              coord1 < halo || coord1 >= size1 + halo ||
                              coord2 < halo || coord2 >= size2 + halo);

        if (!in_halo)
            continue;

        for (int i = 0; i < nfields; i++)
        {
            const long base = fields[i]->lattice().siteFirst() -
                              (long) halo * (fields[i]->lattice().jump(0) + fields[i]->lattice().jump(1) + fields[i]->lattice().jump(2));

            sites[i] = Site(fields[i]->lattice(), base
                            + idx*fields[i]->lattice().jump(0)
                            + coord1*fields[i]->lattice().jump(1)
                            + coord2*fields[i]->lattice().jump(2));
        }

        funct(fields, sites, nfields, params);
    }
}

template <int components = 1>
__host__ __device__ void lattice_add(Field<Real> * fields[], Site * sites, int nfields, double * params, double * outputs)
{
	if constexpr (components == 1)
		(*fields[0])(sites[0]) += (*params);
	else
	{
		#pragma unroll
		for (int i = 0; i < components; i++)
			(*fields[0])(sites[0], i) += (*params);
	}
}

template <int components = 1>
struct lattice_add_functor
{
    __host__ __device__ void operator()(Field<Real> * fields[], Site * sites, int nfields, double * params, double * outputs, void ** vparams = nullptr)
	{
		lattice_add<components>(fields, sites, nfields, params, outputs);
	}
};

template <int components = 1>
__host__ __device__ void lattice_multiply(Field<Real> * fields[], Site * sites, int nfields, double * params, double * outputs)
{
	if constexpr (components == 1)
		(*fields[0])(sites[0]) *= (*params);
	else
	{
		#pragma unroll
		for (int i = 0; i < components; i++)
			(*fields[0])(sites[0], i) *= (*params);
	}
}

template <int components = 1>
struct lattice_multiply_functor
{
    __host__ __device__ void operator()(Field<Real> * fields[], Site * sites, int nfields, double * params, double * outputs, void ** vparams = nullptr)
	{
        (void) vparams;
		lattice_multiply<components>(fields, sites, nfields, params, outputs);
	}
};

#endif // LATTICE_LOOP_HPP