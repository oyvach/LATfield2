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
 		temporaryMemFFT(bool managed = false);
 		~temporaryMemFFT();
 		temporaryMemFFT(long size);

 		int setTemp(long size);
		int reserveDeviceWorkspaceBytes(size_t bytes, const char* context = nullptr, bool allow_shrink = false);
    	void clear();
		long allocated(){return allocated_;}
		size_t deviceAllocated(){return device_allocated_;}
		size_t deviceWorkspaceBytes();
		void * deviceWorkspace();
		size_t minimum_size_;
		size_t current_cap_;
		size_t device_allocated_; //number of complex values in each logical device buffer

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
	// made public for shrinking in main loop
	int reserveDeviceComplexCapacity(size_t capacity, const char* context, bool allow_shrink = false);

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
		bool is_managed;

		size_t deviceComplexBytes();
		void updateDeviceBufferPointers();
		void warnDeviceWorkspaceGrowth(size_t old_bytes, size_t new_bytes, const char* context);
 	};


const int FFT_FORWARD = 1;
const int FFT_BACKWARD = -1;
const int FFT_IN_PLACE = 16;
const int FFT_OUT_OF_PLACE = -16;


#ifndef DOXYGEN_SHOULD_SKIP_THIS

temporaryMemFFT::temporaryMemFFT(bool managed)
{
	temp1_=nullptr;
	temp2_=nullptr;
	device_block_=nullptr;
	temp3_=nullptr;
	temp4_=nullptr;
	temp5_=nullptr;
	allocated_=0;
	device_allocated_=0;
	minimum_size_ = 0;
	is_managed = managed;
	current_cap_ = 0;
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
	minimum_size_ = 0;
	current_cap_ = 0;
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
		          << " bytes). "<< std::endl;
				//   Consider preallocating a larger shared workspace." << std::endl;
	}
}

int temporaryMemFFT::reserveDeviceComplexCapacity(size_t capacity, const char* context, bool allow_shrink)
{
	capacity = std::max(capacity, minimum_size_);
	current_cap_ = std::max(current_cap_, capacity);
	if (capacity <= device_allocated_)
	{
    	if ( 
#ifdef FREETEMPPERCENT
			!allow_shrink || (capacity * 100 > device_allocated_ * (size_t) (100 - FREETEMPPERCENT))
#else
			true
#endif
			)
			return 1;
		else
		{
			// the below will run and reallocate memory to a smaller size
			std::cerr << "Shrinking device workspace capacity from " << device_allocated_ << " to " << capacity
			          << " complex values " << std::endl;
			// size_t minimum_capacity = ...;// previous was wrong – but don't use now anyway
			// capacity = std::max(capacity, minimum_capacity);
		}
	}

	size_t old_bytes = deviceWorkspaceBytes();
	size_t new_bytes = 3 * capacity * deviceComplexBytes();

	if (device_block_ != nullptr) cudaFree(device_block_);

	cudaError_t success;

	if (!is_managed)
	{
		success = cudaMalloc((void **)&device_block_, new_bytes);
	}
	else
	{
		success = cudaMallocManaged((void **)&device_block_, new_bytes);
	}

	if (success != cudaSuccess)
	{
		std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
		std::cerr << "temporaryMemFFT::reserveDeviceComplexCapacity failure details: "
		          << "requested_capacity=" << capacity
		          << ", requested_bytes=" << new_bytes
		          << ", current_cap_=" << current_cap_
		          << ", device_allocated_=" << device_allocated_
		          << ", minimum_size_=" << minimum_size_
		          << ", deviceComplexBytes()=" << deviceComplexBytes()
		          << ", allow_shrink=" << allow_shrink
		          << ", context=" << (context != nullptr ? context : "(null)")
		          << std::endl;
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

int temporaryMemFFT::reserveDeviceWorkspaceBytes(size_t bytes, const char* context, bool allow_shrink)
{
	if (bytes == 0) return 1;

	size_t per_buffer_bytes = 3 * deviceComplexBytes();
	size_t required_capacity = (bytes + per_buffer_bytes - 1) / per_buffer_bytes;
	return reserveDeviceComplexCapacity(required_capacity, context, allow_shrink);
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
		temp1_=nullptr;
		temp2_=nullptr;
		device_block_=nullptr;
		temp3_=nullptr;
		temp4_=nullptr;
		temp5_=nullptr;
}


#endif
//////////////////////Temp memory///////////////////////////

temporaryMemFFT tempMemory(true); // NB: now using managed memory!! experimental




#endif
