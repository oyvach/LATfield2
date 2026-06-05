#ifndef LATFIELD2D_PLANFFT_HPP
#define LATFIELD2D_PLANFFT_HPP

/*! \file LATfield2_PlanFFT.hpp
 \brief FFT wrapper
 LATfield2_PlanFFT.hpp contain the class PlanFFT definition.
 */

//#include "LATfield2_PlanFFT_decl.hpp"

/*! \class temporaryMemFFT
 \brief A class wich handle the additional memory needed by the class PlanFFT_CPU; No documentation!
 */
 class temporaryMemFFT
 	{
 	public:
 		temporaryMemFFT();
 		~temporaryMemFFT();
 		temporaryMemFFT(long size);

		void setDeviceWorkspaceManaged(bool managed);
 		int setTemp(long size);
		int reserveDeviceWorkspaceBytes(size_t bytes, const char* context = nullptr);
    	void clear();
		long allocated(){return allocated_;}
		size_t deviceAllocated(){return device_allocated_;}
		size_t deviceWorkspaceBytes();
		void * deviceWorkspace();

 #ifdef SINGLE
 		fftwf_complex* temp1(){return temp1_;}
 		fftwf_complex* temp2(){return temp2_;}
    cufftComplex * temp3(){return temp3_;}
	cufftComplex * temp4(){return temp4_;}
    cufftComplex * temp5(){return temp5_;}
 #endif

 #ifndef SINGLE
 	        fftw_complex * temp1(){return temp1_;}
 	        fftw_complex * temp2(){return temp2_;}
    cufftDoubleComplex * temp3(){return temp3_;}
	cufftDoubleComplex * temp4(){return temp4_;}
    cufftDoubleComplex * temp5(){return temp5_;}
 #endif

 	private:
#ifdef SINGLE
 		fftwf_complex * temp1_;
 		fftwf_complex * temp2_;
    cufftComplex * device_block_;
    cufftComplex * temp3_;
	cufftComplex * temp4_;
    cufftComplex * temp5_;
 #endif

#ifndef SINGLE
 		fftw_complex * temp1_;
 		fftw_complex * temp2_;
    cufftDoubleComplex * device_block_;
    cufftDoubleComplex * temp3_;
	cufftDoubleComplex * temp4_;
    cufftDoubleComplex * temp5_;
#endif
 		long allocated_; //number of variable stored (bit = allocated*sizeof(fftw(f)_complex))
		size_t device_allocated_; //number of complex values in each logical device buffer
		bool device_workspace_managed_;

		size_t deviceComplexBytes();
		int reserveDeviceComplexCapacity(size_t capacity, const char* context);
		void updateDeviceBufferPointers();
		void warnDeviceWorkspaceGrowth(size_t old_bytes, size_t new_bytes, const char* context);
 	};


const int FFT_FORWARD = 1;
const int FFT_BACKWARD = -1;
const int FFT_IN_PLACE = 16;
const int FFT_OUT_OF_PLACE = -16;


#ifndef DOXYGEN_SHOULD_SKIP_THIS

temporaryMemFFT::temporaryMemFFT()
{
	temp1_=nullptr;
	temp2_=nullptr;
	device_block_=nullptr;
	temp3_=nullptr;
	temp4_=nullptr;
	temp5_=nullptr;
	allocated_=0;
	device_allocated_=0;
	device_workspace_managed_=false;
}
temporaryMemFFT::~temporaryMemFFT()
{
}

temporaryMemFFT::temporaryMemFFT(long size)
{
	temp1_=nullptr;
	temp2_=nullptr;
	device_block_=nullptr;
	temp3_=nullptr;
	temp4_=nullptr;
	temp5_=nullptr;
	allocated_=0;
	device_allocated_=0;
	device_workspace_managed_=false;
	setTemp(size);
}

size_t temporaryMemFFT::deviceComplexBytes()
{
#ifdef SINGLE
	return sizeof(cufftComplex);
#else
	return sizeof(cufftDoubleComplex);
#endif
}

size_t temporaryMemFFT::deviceWorkspaceBytes()
{
	return 3 * device_allocated_ * deviceComplexBytes();
}

void * temporaryMemFFT::deviceWorkspace()
{
	return reinterpret_cast<void*>(device_block_);
}

void temporaryMemFFT::updateDeviceBufferPointers()
{
	if (device_block_ == nullptr)
	{
		temp3_ = nullptr;
		temp4_ = nullptr;
		temp5_ = nullptr;
		return;
	}

	temp3_ = device_block_;
	temp4_ = device_block_ + device_allocated_;
	temp5_ = device_block_ + 2 * device_allocated_;
}

void temporaryMemFFT::warnDeviceWorkspaceGrowth(size_t old_bytes, size_t new_bytes, const char* context)
{
	if (old_bytes == 0 || context == nullptr) return;

	int rank = 0;
#ifdef MPI_VERSION
	int mpi_initialized = 0;
	MPI_Initialized(&mpi_initialized);
	if (mpi_initialized) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	if (rank == 0)
	{
		std::cerr << "LATfield2 temporary device workspace growth during " << context
		          << " (old=" << old_bytes << " bytes, new=" << new_bytes
		          << " bytes). Consider preallocating a larger shared workspace." << std::endl;
	}
}

void temporaryMemFFT::setDeviceWorkspaceManaged(bool managed)
{
	if (!managed || device_workspace_managed_) return;

	size_t old_capacity = device_allocated_;
	device_workspace_managed_ = true;

	if (device_block_ == nullptr) return;

	cudaFree(device_block_);
	device_block_ = nullptr;
	device_allocated_ = 0;
	updateDeviceBufferPointers();
	reserveDeviceComplexCapacity(old_capacity, "temporaryMemFFT::setDeviceWorkspaceManaged");
}

int temporaryMemFFT::reserveDeviceComplexCapacity(size_t capacity, const char* context)
{
	if (capacity <= device_allocated_) return 1;

	size_t old_bytes = deviceWorkspaceBytes();
	size_t new_bytes = 3 * capacity * deviceComplexBytes();

	if (device_block_ != nullptr) cudaFree(device_block_);

	auto success = device_workspace_managed_
		? cudaMallocManaged((void **)&device_block_, new_bytes)
		: cudaMalloc((void **)&device_block_, new_bytes);

	if (success != cudaSuccess)
	{
		std::cerr << (device_workspace_managed_ ? "cudaMallocManaged" : "cudaMalloc")
		          << " failed: " << cudaGetErrorString(success) << std::endl;
		device_block_ = nullptr;
		device_allocated_ = 0;
		updateDeviceBufferPointers();
		throw std::runtime_error("Memory allocation failed in temporaryMemFFT::reserveDeviceComplexCapacity");
	}

	device_allocated_ = capacity;
	updateDeviceBufferPointers();
	warnDeviceWorkspaceGrowth(old_bytes, new_bytes, context);
	return 1;
}

int temporaryMemFFT::reserveDeviceWorkspaceBytes(size_t bytes, const char* context)
{
	if (bytes == 0) return 1;

	size_t per_buffer_bytes = 3 * deviceComplexBytes();
	size_t required_capacity = (bytes + per_buffer_bytes - 1) / per_buffer_bytes;
	return reserveDeviceComplexCapacity(required_capacity, context);
}

int temporaryMemFFT::setTemp(long size)
{
	if(size>allocated_)
	{
#ifdef SINGLE
		if(allocated_)
		{
			if(temp1_!=nullptr)fftwf_free(temp1_);
			if(temp2_!=nullptr)fftwf_free(temp2_);
		}
		temp1_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
		temp2_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));

		//debug cout<<"("<< parallel.grid_rank()[0]<< ","<<parallel.grid_rank()[1] <<"): called temporary resize. Old size: " <<allocated_<<" , new size: "<< size<<endl;

#endif

#ifndef SINGLE
		if(allocated_)
		{
			if(temp1_!=nullptr)fftw_free(temp1_);
			if(temp2_!=nullptr)fftw_free(temp2_);
		}
		temp1_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
		temp2_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
#endif

		if (temp1_ == nullptr || temp2_ == nullptr)
		{
			throw std::runtime_error("Host memory allocation failed in temporaryMemFFT::setTemp");
		}

		allocated_ = size ;
	}
	reserveDeviceComplexCapacity((size_t)size, nullptr);
	return 1;
}

void temporaryMemFFT::clear()
{
#ifdef SINGLE
		if(temp1_!=nullptr)fftwf_free(temp1_);
		if(temp2_!=nullptr)fftwf_free(temp2_);
#endif
#ifndef SINGLE
		if(temp1_!=nullptr)fftw_free(temp1_);
		if(temp2_!=nullptr)fftw_free(temp2_);
#endif
		if(device_block_!=nullptr)cudaFree(device_block_);
		allocated_ = 0;
		device_allocated_ = 0;
		device_workspace_managed_ = false;
		temp1_=nullptr;
		temp2_=nullptr;
		device_block_=nullptr;
		temp3_=nullptr;
		temp4_=nullptr;
		temp5_=nullptr;
}


#endif
//////////////////////Temp memory///////////////////////////

temporaryMemFFT tempMemory;




#endif
