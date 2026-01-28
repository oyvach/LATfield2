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

 		int setTemp(long size);
    	void clear();
		long allocated(){return allocated_;}

 #ifdef SINGLE
 		fftwf_complex* temp1(){return temp1_;}
 		fftwf_complex* temp2(){return temp2_;}
    cufftComplex * temp3(){return temp3_;}
	cufftComplex * temp4(){return temp4_;}
 #endif

 #ifndef SINGLE
 	        fftw_complex * temp1(){return temp1_;}
 	        fftw_complex * temp2(){return temp2_;}
    cufftDoubleComplex * temp3(){return temp3_;}
	cufftDoubleComplex * temp4(){return temp4_;}
 #endif

 	private:
 #ifdef SINGLE
 		fftwf_complex * temp1_;
 		fftwf_complex * temp2_;
    cufftComplex * temp3_;
	cufftComplex * temp4_;
 #endif

 #ifndef SINGLE
 		fftw_complex * temp1_;
 		fftw_complex * temp2_;
    cufftDoubleComplex * temp3_;
	cufftDoubleComplex * temp4_;
 #endif
 		long allocated_; //number of variable stored (bit = allocated*sizeof(fftw(f)_complex))
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
	temp3_=nullptr;
	temp4_=nullptr;
	allocated_=0;
}
temporaryMemFFT::~temporaryMemFFT()
{
}

temporaryMemFFT::temporaryMemFFT(long size)
{
#ifdef SINGLE
	temp1_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
	temp2_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
	auto success = cudaMalloc(&temp3_, size * sizeof(cufftComplex));
#endif

#ifndef SINGLE
	temp1_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	temp2_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	auto success = cudaMalloc(&temp3_, size * sizeof(cufftDoubleComplex));
#endif

	if (success != cudaSuccess)
    {
        std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
        throw std::runtime_error("Memory allocation failed in constructor of temporaryMemFFT");
    }

#ifdef SINGLE
	success = cudaMalloc(&temp4_, size * sizeof(cufftComplex));
#else
	success = cudaMalloc(&temp4_, size * sizeof(cufftDoubleComplex));
#endif

	if (success != cudaSuccess)
	{
		std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
		throw std::runtime_error("Memory allocation failed in constructor of temporaryMemFFT");
	}
	allocated_=size;
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
			if(temp3_!=nullptr)cudaFree(temp3_);
			if(temp4_!=nullptr)cudaFree(temp4_);
		}
		temp1_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
		temp2_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
		auto success = cudaMalloc(&temp3_, size * sizeof(cufftComplex));

		//debug cout<<"("<< parallel.grid_rank()[0]<< ","<<parallel.grid_rank()[1] <<"): called temporary resize. Old size: " <<allocated_<<" , new size: "<< size<<endl;

#endif

#ifndef SINGLE
		if(allocated_)
		{
			if(temp1_!=nullptr)fftw_free(temp1_);
			if(temp2_!=nullptr)fftw_free(temp2_);
			if(temp3_!=nullptr)cudaFree(temp3_);
			if(temp4_!=nullptr)cudaFree(temp4_);
		}
		temp1_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
		temp2_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
		auto success = cudaMalloc(&temp3_, size * sizeof(cufftDoubleComplex));
#endif

		if (success != cudaSuccess)
		{
			std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
			throw std::runtime_error("Memory allocation failed in temporaryMemFFT::setTemp");
		}

#ifdef SINGLE
		success = cudaMalloc(&temp4_, size * sizeof(cufftComplex));
#else
		success = cudaMalloc(&temp4_, size * sizeof(cufftDoubleComplex));
#endif

		if (success != cudaSuccess)
		{
			std::cerr << "cudaMalloc failed: " << cudaGetErrorString(success) << std::endl;
			throw std::runtime_error("Memory allocation failed in temporaryMemFFT::setTemp");
		}

		allocated_ = size ;
	}
	return 1;
}

void temporaryMemFFT::clear()
{
	if(allocated_>0)
	{
#ifdef SINGLE
		if(temp1_!=nullptr)fftwf_free(temp1_);
		if(temp2_!=nullptr)fftwf_free(temp2_);
#endif
#ifndef SINGLE
		if(temp1_!=nullptr)fftw_free(temp1_);
		if(temp2_!=nullptr)fftw_free(temp2_);
#endif
		if(temp3_!=nullptr)cudaFree(temp3_);
		if(temp4_!=nullptr)cudaFree(temp4_);
		allocated_ = 0;
		temp1_=nullptr;
		temp2_=nullptr;
		temp3_=nullptr;
		temp4_=nullptr;
	}
}


#endif
//////////////////////Temp memory///////////////////////////

temporaryMemFFT tempMemory;




#endif
