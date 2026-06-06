#ifndef LATFIELD2_PLANFFT_DECL_HPP
#define LATFIELD2_PLANFFT_DECL_HPP

#include <cstdlib>
#include <cstring>

#ifdef SINGLE
#define MPI_DATA_PREC MPI_FLOAT
#endif

#ifndef SINGLE
#define MPI_DATA_PREC MPI_DOUBLE
#endif

#ifndef NULLFFTWPLAN
#ifndef SINGLE
#define NULLFFTWPLAN static_cast<fftw_plan>(NULL)
#else
#define NULLFFTWPLAN static_cast<fftwf_plan>(NULL)
#endif

#endif


/* these are not the actual compiler variables, these are
   only the declaration that they exist somewhere at all. */

extern  const int FFT_FORWARD;
extern  const int FFT_BACKWARD;
extern  const int FFT_IN_PLACE;
extern  const int FFT_OUT_OF_PLACE;
extern  const int FFT_EXECUTION_HOST_MPI;
extern  const int FFT_EXECUTION_CUDA_AWARE_MPI;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

const int FFT_EXECUTION_HOST_MPI = 0;
const int FFT_EXECUTION_CUDA_AWARE_MPI = 1;

const char* cufftGetErrorString(cufftResult status) {
    switch (status) {
        case CUFFT_SUCCESS:                return "CUFFT_SUCCESS";
        case CUFFT_INVALID_PLAN:           return "CUFFT_INVALID_PLAN";
        case CUFFT_ALLOC_FAILED:           return "CUFFT_ALLOC_FAILED";
        case CUFFT_INVALID_TYPE:           return "CUFFT_INVALID_TYPE";
        case CUFFT_INVALID_VALUE:          return "CUFFT_INVALID_VALUE";
        case CUFFT_INTERNAL_ERROR:         return "CUFFT_INTERNAL_ERROR";
        case CUFFT_EXEC_FAILED:            return "CUFFT_EXEC_FAILED";
        case CUFFT_SETUP_FAILED:           return "CUFFT_SETUP_FAILED";
        case CUFFT_INVALID_SIZE:           return "CUFFT_INVALID_SIZE";
        case CUFFT_UNALIGNED_DATA:         return "CUFFT_UNALIGNED_DATA";
#if (CUDART_VERSION >= 6050)
        case CUFFT_INCOMPLETE_PARAMETER_LIST: return "CUFFT_INCOMPLETE_PARAMETER_LIST";
        case CUFFT_INVALID_DEVICE:         return "CUFFT_INVALID_DEVICE";
        case CUFFT_PARSE_ERROR:            return "CUFFT_PARSE_ERROR";
        case CUFFT_NO_WORKSPACE:           return "CUFFT_NO_WORKSPACE";
        case CUFFT_NOT_IMPLEMENTED:        return "CUFFT_NOT_IMPLEMENTED";
        case CUFFT_LICENSE_ERROR:          return "CUFFT_LICENSE_ERROR";
#endif
        default:                           return "Unknown CUFFT error";
    }
}

extern  temporaryMemFFT tempMemory;

/*! \class PlanFFT

 \brief Class which handle fourier transforms of fields on 3d cubic lattices.
 This class allow to perform fourier transform of real and complex fields. See poissonSolver example to have have a short intro of usage.

 One should understand that first a plan is created then execute (in the FFTW fashion). The plan link to fields, one on fourier space, one on real space. Both field will be allocated by the planer. But need to be initialized.

 One need to be carefull to corretly define the lattice and field.
 \sa void Lattice::initializeRealFFT(Lattice & lat_real, int halo);
 \sa void Lattice::initializeComplexFFT(Lattice & lat_real, int halo);
 For more detail see the QuickStart guide.
 */
template<class compType>
class PlanFFT
{
public:
  //! Constructor.
  PlanFFT();

  //! Destructor.
  ~PlanFFT();

#ifndef SINGLE
  /*!
   Constructor with initialization for complex to complex tranform.
   \sa initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);
   \param rfield : real space field
   \param kfield : fourier space field
   \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both fourier and real space field point to the same data array.
   */
  PlanFFT(Field<compType>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE, const bool managed = false);
  /*!
   initialization for complex to complex tranform.
   For more detail see the QuickStart guide.
   \param rfield : real space field
   \param kfield : fourier space field
   \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both fourier and real space field point to the same data array.
   */
  void initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);

  /*!
   Constructor with initialization for real to complex tranform.
   \sa initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);
   \param rfield : real space field
   \param kfield : fourier space field
   \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both fourier and real space field point to the same data array.
   */
  PlanFFT(Field<double>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE, const bool managed = false);
  /*!
   initialization for real to complex tranform.
   For more detail see the QuickStart guide.
   \param rfield : real space field
   \param kfield : fourier space field
   \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both fourier and real space field point to the same data array.
   */
  void initialize(Field<double>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);



#endif

#ifdef SINGLE

  PlanFFT(Field<compType>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE, const bool managed = false);
  void initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);


  PlanFFT(Field<float>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE, const bool managed = false);
  void initialize(Field<float>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);


#endif


  void preallocate();
  void preallocate(size_t min_shared_device_workspace_bytes);
  void execute(int fft_type);
  void setExecutionMode(int mode);
  int executionMode() const;
  bool usingCudaAwareMPI() const;

private:
  void PrintPlans() {
#ifndef SINGLE
    std::cout << fPlan_i_ << " "; fftw_print_plan(fPlan_i_); std::cout << std::endl;
    std::cout << fPlan_j_ << " "; fftw_print_plan(fPlan_j_); std::cout << std::endl;
    std::cout << fPlan_k_ << " "; fftw_print_plan(fPlan_k_); std::cout << std::endl;
    std::cout << fPlan_k_real_ << " "; fftw_print_plan(fPlan_k_real_); std::cout << std::endl;
    std::cout << bPlan_k_ << " "; fftw_print_plan(bPlan_k_); std::cout << std::endl;
    std::cout << bPlan_j_ << " "; fftw_print_plan(bPlan_j_); std::cout << std::endl;
    std::cout << bPlan_j_real_ << " "; fftw_print_plan(bPlan_j_real_); std::cout << std::endl;
    std::cout << bPlan_i_ << " "; fftw_print_plan(bPlan_i_); std::cout << std::endl;
#else
    std::cout << fPlan_i_ << " "; fftwf_print_plan(fPlan_i_); std::cout << std::endl;
    std::cout << fPlan_j_ << " "; fftwf_print_plan(fPlan_j_); std::cout << std::endl;
    std::cout << fPlan_k_ << " "; fftwf_print_plan(fPlan_k_); std::cout << std::endl;
    std::cout << fPlan_k_real_ << " "; fftwf_print_plan(fPlan_k_real_); std::cout << std::endl;
    std::cout << bPlan_k_ << " "; fftwf_print_plan(bPlan_k_); std::cout << std::endl;
    std::cout << bPlan_j_ << " "; fftwf_print_plan(bPlan_j_); std::cout << std::endl;
    std::cout << bPlan_j_real_ << " "; fftwf_print_plan(bPlan_j_real_); std::cout << std::endl;
    std::cout << bPlan_i_ << " "; fftwf_print_plan(bPlan_i_); std::cout << std::endl;
#endif


}

private:
  bool status_;
  bool type_;
  int mem_type_;
  bool managed_;
  bool alignment_even_;
  int execution_mode_;
  bool cuda_aware_mpi_enabled_;
  bool cuda_aware_mpi_resolved_;
  bool cuda_aware_mpi_warning_emitted_;

  static bool R2C;
  static bool C2C;

  //static bool initialized;

  //data description variable, fftw plan, and temp

  int components_;
  int rSize_[3];
  int kSize_[3];
  int rJump_[3];
  int kJump_[3];
  int rSizeLocal_[3];
  int kSizeLocal_[3];
  int r2cSize_;  //of 0 dimension
  int r2cSizeLocal_; //of r2csize on proc_dim[1]
  int r2cSizeLocal_as_;
  int rHalo_;
  int kHalo_;

  void resolveCudaAwareMPI();
  bool envVarEnabled(const char* name) const;
  void warnCudaAwareMPIFallback(const char* reason);
  long requiredTemporaryCapacity() const;

  cufftHandle cufPlan_i_;
  cufftHandle cufPlan_j_;
  cufftHandle cufPlan_k_;
  cufftHandle cufPlan_k_real_;

  cufftHandle cubPlan_i_;
  cufftHandle cubPlan_j_;
  cufftHandle cubPlan_j_real_;
  cufftHandle cubPlan_k_;

#ifdef SINGLE
  float * rData_; //pointer to start of data (halo skip)
  fftwf_complex * cData_; //pointer to start of data (halo skip)
  fftwf_complex * kData_; //pointer to start of data (halo skip)
  fftwf_complex * temp_;
  fftwf_complex * temp1_;//needed if field got more than 1 component

  fftwf_plan fPlan_i_;
  fftwf_plan fPlan_j_;
  fftwf_plan fPlan_k_;
  fftwf_plan fPlan_k_real_;

  fftwf_plan bPlan_i_;
  fftwf_plan bPlan_j_;
  fftwf_plan bPlan_j_real_;
  fftwf_plan bPlan_k_;

  ///transpostion fonction

  /// forward real to complex
  // first transopsition
  void transpose_0_2( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k);
  void transpose_0_2( cufftComplex * in, cufftComplex * out,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream);
  void transpose_0_2_last_proc( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k);
  void transpose_0_2_last_proc( cufftComplex * in, cufftComplex * out,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream);
  void implement_local_0_last_proc( fftwf_complex * in, fftwf_complex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size);
  void implement_local_0_last_proc( cufftComplex * in, cufftComplex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size, cudaStream_t &stream);
  // second transposition
  void transpose_1_2(fftwf_complex * in , fftwf_complex * out  ,int dim_i,int dim_j ,int dim_k);
  void transpose_1_2(cufftComplex * in , cufftComplex * out  ,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream);
  //third transposition
  void transpose_back_0_3(fftwf_complex * in, fftwf_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components,int comp);
  void implement_0(fftwf_complex * in, fftwf_complex * out,int r2c_size,int local_size_j,int local_size_k,int halo,int components,int comp);
  ////backward real to complex
  void b_arrange_data_0(fftwf_complex *in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k, int khalo, int components, int comp);
  void b_transpose_back_0_1(fftwf_complex * in, fftwf_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size);
  void b_transpose_back_0_1(cufftComplex * in, cufftComplex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size, cudaStream_t &stream);
  void b_implement_0(fftwf_complex * in, fftwf_complex * out,int r2c_size,int local_size_j,int local_size_k);
  void b_implement_0(cufftComplex * in, cufftComplex * out,int r2c_size,int local_size_j,int local_size_k, cudaStream_t &stream);

#endif
#ifndef SINGLE


  double * rData_; //pointer to start of data (halo skip)
  fftw_complex * cData_; //pointer to start of data (halo skip)
  fftw_complex * kData_; //pointer to start of data (halo skip)
  fftw_complex * temp_;
  fftw_complex * temp1_;//needed if field got more than 1 component


  fftw_plan fPlan_i_;
  fftw_plan fPlan_j_;
  fftw_plan fPlan_k_;
  fftw_plan fPlan_k_real_;

  fftw_plan bPlan_i_;
  fftw_plan bPlan_j_;
  fftw_plan bPlan_j_real_;
  fftw_plan bPlan_k_;

  ///transpostion fonction

  /// forward real to complex
  // first transopsition
  void transpose_0_2( fftw_complex * in, fftw_complex * out,int dim_i,int dim_j ,int dim_k);
  void transpose_0_2( cufftDoubleComplex * in, cufftDoubleComplex * out,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream);
  void transpose_0_2_last_proc( fftw_complex * in, fftw_complex * out,int dim_i,int dim_j ,int dim_k);
  void transpose_0_2_last_proc( cufftDoubleComplex * in, cufftDoubleComplex * out,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream);
  void implement_local_0_last_proc( fftw_complex * in, fftw_complex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size);
  void implement_local_0_last_proc( cufftDoubleComplex * in, cufftDoubleComplex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size, cudaStream_t &stream);
  // second transposition
  void transpose_1_2(fftw_complex * in , fftw_complex * out  ,int dim_i,int dim_j ,int dim_k);
  void transpose_1_2(cufftDoubleComplex * in , cufftDoubleComplex * out  ,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream);
  //third transposition
  void transpose_back_0_3(fftw_complex * in, fftw_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components,int comp);
  void implement_0(fftw_complex * in, fftw_complex * out,int r2c_size,int local_size_j,int local_size_k,int halo,int components,int comp);
  ////backward real to complex
  void b_arrange_data_0(fftw_complex *in, fftw_complex * out,int dim_i,int dim_j ,int dim_k, int khalo, int components, int comp);
  void b_transpose_back_0_1(fftw_complex * in, fftw_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size);
  void b_transpose_back_0_1(cufftDoubleComplex * in, cufftDoubleComplex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size, cudaStream_t &stream);
  void b_implement_0(fftw_complex * in, fftw_complex * out,int r2c_size,int local_size_j,int local_size_k);
  void b_implement_0(cufftDoubleComplex * in, cufftDoubleComplex * out,int r2c_size,int local_size_j,int local_size_k, cudaStream_t &stream);




#endif

};

//constants
//template<class compType>
//bool PlanFFT<compType>::initialized = true;
template<class compType>
bool PlanFFT<compType>::R2C=false;
template<class compType>
bool PlanFFT<compType>::C2C=true;



template<class compType>
PlanFFT<compType>::~PlanFFT() {
#ifndef SINGLE
  //if (fPlan_i_ != NULLFFTWPLAN) { fftw_destroy_plan(fPlan_i_); }
  //if (fPlan_j_ != NULLFFTWPLAN) { fftw_destroy_plan(fPlan_j_); }
  //if (fPlan_k_ != NULLFFTWPLAN) { fftw_destroy_plan(fPlan_k_); }
  //if (fPlan_k_real_ != NULLFFTWPLAN) { fftw_destroy_plan(fPlan_k_real_); }
  //if (bPlan_k_ != NULLFFTWPLAN) { fftw_destroy_plan(bPlan_k_); }
  //if (bPlan_j_ != NULLFFTWPLAN) { fftw_destroy_plan(bPlan_j_); }
  //if (bPlan_j_real_ != NULLFFTWPLAN) { fftw_destroy_plan(bPlan_j_real_); }
  //if (bPlan_i_ != NULLFFTWPLAN) { fftw_destroy_plan(bPlan_i_); }
#else
  //if (fPlan_i_ != NULLFFTWPLAN) { fftwf_destroy_plan(fPlan_i_); }
  //if (fPlan_j_ != NULLFFTWPLAN) { fftwf_destroy_plan(fPlan_j_); }
  //if (fPlan_k_ != NULLFFTWPLAN) { fftwf_destroy_plan(fPlan_k_); }
  //if (fPlan_k_real_ != NULLFFTWPLAN) { fftwf_destroy_plan(fPlan_k_real_); }
  //if (bPlan_k_ != NULLFFTWPLAN) { fftwf_destroy_plan(bPlan_k_); }
  //if (bPlan_j_ != NULLFFTWPLAN) { fftwf_destroy_plan(bPlan_j_); }
  //if (bPlan_j_real_ != NULLFFTWPLAN) { fftwf_destroy_plan(bPlan_j_real_); }
  //if (bPlan_i_ != NULLFFTWPLAN) { fftwf_destroy_plan(bPlan_i_); }
#endif



}


template<class compType>
PlanFFT<compType>::PlanFFT() :
fPlan_i_(NULLFFTWPLAN),
fPlan_j_(NULLFFTWPLAN),
fPlan_k_(NULLFFTWPLAN),
fPlan_k_real_(NULLFFTWPLAN),
bPlan_k_(NULLFFTWPLAN),
bPlan_j_(NULLFFTWPLAN),
bPlan_j_real_(NULLFFTWPLAN),
bPlan_i_(NULLFFTWPLAN)
{
  status_ = false;
  alignment_even_ = true;
  execution_mode_ = FFT_EXECUTION_HOST_MPI;
  cuda_aware_mpi_enabled_ = false;
  cuda_aware_mpi_resolved_ = false;
  cuda_aware_mpi_warning_emitted_ = false;
}

template<class compType>
void PlanFFT<compType>::setExecutionMode(int mode)
{
  execution_mode_ = mode;
  cuda_aware_mpi_enabled_ = false;
  cuda_aware_mpi_resolved_ = false;
}

template<class compType>
int PlanFFT<compType>::executionMode() const
{
  return execution_mode_;
}

template<class compType>
bool PlanFFT<compType>::usingCudaAwareMPI() const
{
  return cuda_aware_mpi_enabled_;
}

template<class compType>
bool PlanFFT<compType>::envVarEnabled(const char* name) const
{
  const char* value = std::getenv(name);
  if (value == nullptr) return false;
  return !(strcmp(value, "0") == 0 || strcmp(value, "false") == 0 || strcmp(value, "FALSE") == 0);
}

template<class compType>
void PlanFFT<compType>::warnCudaAwareMPIFallback(const char* reason)
{
  if (cuda_aware_mpi_warning_emitted_) return;
  if (parallel.isRoot())
  {
    COUT << "PlanFFT: disabling experimental CUDA-aware MPI path, falling back to host-staged MPI";
    if (reason != nullptr) COUT << " (" << reason << ")";
    COUT << endl;
  }
  cuda_aware_mpi_warning_emitted_ = true;
}

template<class compType>
void PlanFFT<compType>::resolveCudaAwareMPI()
{
  if (cuda_aware_mpi_resolved_) return;

  cuda_aware_mpi_resolved_ = true;
  cuda_aware_mpi_enabled_ = false;

  if (execution_mode_ != FFT_EXECUTION_CUDA_AWARE_MPI) return;

  if (envVarEnabled("LATFIELD2_DISABLE_CUDA_AWARE_MPI"))
  {
    warnCudaAwareMPIFallback("disabled by LATFIELD2_DISABLE_CUDA_AWARE_MPI");
    return;
  }

  if (tempMemory.temp5() == nullptr)
  {
    warnCudaAwareMPIFallback("device communication buffer not allocated");
    return;
  }

  const bool explicit_enable =
    envVarEnabled("LATFIELD2_ENABLE_CUDA_AWARE_MPI") ||
    envVarEnabled("MPICH_GPU_SUPPORT_ENABLED") ||
    envVarEnabled("MV2_USE_CUDA") ||
    envVarEnabled("PSM2_CUDA") ||
    envVarEnabled("OMPI_MCA_opal_cuda_support") ||
    envVarEnabled("OMPI_MCA_mpi_cuda_support");

  if (!explicit_enable)
  {
    warnCudaAwareMPIFallback("no CUDA-aware MPI support hint found");
    return;
  }

  cuda_aware_mpi_enabled_ = true;
}


#ifdef SINGLE



template<class compType>
PlanFFT<compType>::PlanFFT(Field<compType>* rfield, Field<compType>* kfield,const int mem_type , const bool managed) : PlanFFT()
{
  managed_ = managed;
  status_ = false;
  initialize(rfield, kfield, mem_type);
}

template<class compType>
void PlanFFT<compType>::initialize(Field<compType>*  rfield,Field<compType>*  kfield,const int mem_type )
{
  type_ = C2C;
  mem_type_=mem_type;

  //general variable

  COUT<<"INITIALIZING COMPLEX FFT"<<endl;

  if(rfield->components() != kfield->components())
  {
    cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for fields with same number of components"<<endl;
    cerr<<"Latfield2d::PlanFFT::initialize : coordinate and fourier space fields have not the same number of components"<<endl;
    cerr<<"Latfield2d : Abort Process Requested"<<endl;

  }
  else components_ = rfield->components();

  for(int i = 0; i<3; i++)
  {
    rSize_[i]=rfield->lattice().size(i);
    kSize_[i]=kfield->lattice().size(i);
    rSizeLocal_[i]=rfield->lattice().sizeLocal(i);
    kSizeLocal_[i]=kfield->lattice().sizeLocal(i);
    rJump_[i]=rfield->lattice().jump(i);
    kJump_[i]=kfield->lattice().jump(i);
  }
  r2cSize_ = 0;
  r2cSizeLocal_as_ = 0;
  r2cSizeLocal_ = 0;
  rHalo_ = rfield->lattice().halo();
  kHalo_ = kfield->lattice().halo();

  /////from latfield2d_IO
  tempMemory.setTemp((long)(rSize_[0] + 2*rHalo_)  * (long)(rSizeLocal_[1] + 2*rHalo_) * (long)(rSizeLocal_[2] + 2*rHalo_));

  	temp_  = tempMemory.temp1();
  	temp1_ = tempMemory.temp2();

  	if(rfield->lattice().dim()!=3)
  	{
  		if(parallel.isRoot())
  		{
  			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
  			cerr<<"Latfield2d::PlanFFT::initialize : real lattice have not 3 dimensions"<<endl;
  			cerr<<"Latfield2d : Abort Process Requested"<<endl;

  		}
  		parallel.abortForce();
  	}


  	if(rSize_[0]!=rSize_[1] | rSize_[1]!=rSize_[2])
  	{
  		if(parallel.isRoot())
  		{
  			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
  			cerr<<"Latfield2d::PlanFFT::initialize : real lattice is not cubic"<<endl;
  			cerr<<"Latfield2d : Abort Process Requested"<<endl;

   		}
  		parallel.abortForce();
  	}

  	//initialization of fftw plan

    //Pointer to data

    long rfield_size = rfield->lattice().sitesLocalGross();
  	long kfield_size = kfield->lattice().sitesLocalGross();

  	if(mem_type_ == FFT_IN_PLACE)
  	{
  		if(rfield_size>=kfield_size)
  		{
  			rfield->alloc();
  			kfield->data() = (Imag*)rfield->data();
  		}
  		else
  		{
  			kfield->alloc();
  			rfield->data() = (Imag*)kfield->data();
  		}
  	}
  	if(mem_type_ == FFT_OUT_OF_PLACE)
  	{
  		rfield->alloc();
  		kfield->alloc();
  	}

  	rData_ = (float*)rfield->data(); //to be sure that rData is instantiate !
  	cData_ = (fftwf_complex*)rfield->data();
  	cData_ += rfield->lattice().siteFirst()*components_;
  	kData_ = (fftwf_complex*)kfield->data();
  	kData_ += kfield->lattice().siteFirst()*components_;

  	//Forward plan
  	fPlan_i_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,cData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  	fPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE);
  	fPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,kData_,NULL,components_, rJump_[1]*components_,FFTW_FORWARD,FFTW_ESTIMATE);
  	//Backward plan

  	bPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,kData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  	bPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE);
  	bPlan_i_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,cData_,NULL,components_, rJump_[1]*components_,FFTW_BACKWARD,FFTW_ESTIMATE| FFTW_PRESERVE_INPUT);

  	//allocation of field






  ///end of from latfield2d_IO



}

template<class compType>
PlanFFT<compType>::PlanFFT(Field<float>* rfield, Field<compType>*  kfield,const int mem_type, const bool managed )  : PlanFFT()
{
  managed_ = managed;
  status_ = false;
  initialize(rfield,kfield,mem_type);
}

template<class compType>
void PlanFFT<compType>::initialize(Field<float>*  rfield,Field<compType>*   kfield, const int mem_type )
{
  type_ = R2C;
  mem_type_=mem_type;

  //general variable
  if(rfield->components() != kfield->components())
  {
    cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for fields with same number of components"<<endl;
    cerr<<"Latfield2d::PlanFFT::initialize : coordinate and fourier space fields have not the same number of components"<<endl;
    cerr<<"Latfield2d : Abort Process Requested"<<endl;

  }
  else components_ = rfield->components();

  for(int i = 0; i<3; i++)
  {
    rSize_[i]=rfield->lattice().size(i);
    kSize_[i]=kfield->lattice().size(i);
    rSizeLocal_[i]=rfield->lattice().sizeLocal(i);
    kSizeLocal_[i]=kfield->lattice().sizeLocal(i);
    rJump_[i]=rfield->lattice().jump(i);
    kJump_[i]=kfield->lattice().jump(i);
  }
  r2cSize_ = rfield->lattice().size(0)/2 + 1;
  r2cSizeLocal_as_ = rfield->lattice().sizeLocal(0)/(2*parallel.grid_size()[1]);
  if(parallel.last_proc()[1]) r2cSizeLocal_ = r2cSizeLocal_as_ + 1;
  else r2cSizeLocal_ = r2cSizeLocal_as_;
  rHalo_ = rfield->lattice().halo();
  kHalo_ = kfield->lattice().halo();





  preallocate();



  temp_  = tempMemory.temp1();
  temp1_ = tempMemory.temp2();




  if(rfield->lattice().dim()!=3)
  {
    if(parallel.isRoot())
    {
      cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
      cerr<<"Latfield2d::PlanFFT::initialize : real lattice have not 3 dimensions"<<endl;
      cerr<<"Latfield2d : Abort Process Requested"<<endl;

    }
    parallel.abortForce();
  }

  //COUT << rSize_[0] << "   "<< rSize_[1]<< "   "<< rSize_[2] <<endl;

  if(rSize_[0]!=rSize_[1] | rSize_[1]!=rSize_[2])
  {
    if(parallel.isRoot())
    {
      cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
      cerr<<"Latfield2d::PlanFFT::initialize : real lattice is not cubic"<<endl;
      cerr<<"Latfield2d : Abort Process Requested"<<endl;

    }
    parallel.abortForce();
  }

  //create the fftw_plan

  cufftResult cufft_status;
  int inembed[1] = {rSizeLocal_[0]};
  int onembed[1] = {(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2] + 1};

  cufft_status = cufftPlanMany(&cufPlan_i_, 1, &rSize_[0], inembed, components_, rJump_[1]*components_, onembed, rSizeLocal_[1]*rSizeLocal_[2], 1, CUFFT_R2C, rSizeLocal_[1]);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for forward plan i" << endl;
    parallel.abortForce();
  }

  inembed[0] = 1;
  onembed[0] = 1;

  cufft_status = cufftPlanMany(&cufPlan_j_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, 1, onembed, rSizeLocal_[2]*r2cSizeLocal_, 1, CUFFT_C2C, rSizeLocal_[2]*r2cSizeLocal_);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for forward plan j" << endl;
    parallel.abortForce();
  }

  cufft_status = cufftPlanMany(&cufPlan_k_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, 1, onembed, rSizeLocal_[2]*r2cSizeLocal_as_, 1, CUFFT_C2C, r2cSizeLocal_as_);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for forward plan k" << endl;
    parallel.abortForce();
  }

  inembed[0] = r2cSizeLocal_;

  cufft_status = cufftPlanMany(&cufPlan_k_real_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, r2cSizeLocal_, onembed, rSizeLocal_[2], 1, CUFFT_C2C, rSizeLocal_[2]);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for forward plan k_real" << endl;
    parallel.abortForce();
  }

  inembed[0] = 1;

  cufft_status = cufftPlanMany(&cubPlan_k_, 1, &rSize_[0], inembed, kSizeLocal_[2]*r2cSizeLocal_, 1, onembed, kSizeLocal_[2]*r2cSizeLocal_, 1, CUFFT_C2C, kSizeLocal_[2]*r2cSizeLocal_);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for backward plan k" << endl;
    parallel.abortForce();
  }

  cufft_status = cufftPlanMany(&cubPlan_j_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, 1, onembed, rSizeLocal_[2]*r2cSizeLocal_as_, 1, CUFFT_C2C, r2cSizeLocal_as_);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for backward plan j" << endl;
    parallel.abortForce();
  }

  inembed[0] = r2cSizeLocal_;

  cufft_status = cufftPlanMany(&cubPlan_j_real_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, r2cSizeLocal_, onembed, rSizeLocal_[2], 1, CUFFT_C2C, rSizeLocal_[2]);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for backward plan j_real" << endl;
    parallel.abortForce();
  }

  inembed[0] = r2cSize_;
  onembed[0] = rJump_[1];

  cufft_status = cufftPlanMany(&cubPlan_i_, 1, &rSize_[0], inembed, 1, r2cSize_, onembed, components_, rJump_[1]*components_, CUFFT_C2R, rSizeLocal_[1]);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for backward plan i" << endl;
    parallel.abortForce();
  }

  fPlan_i_ = fftwf_plan_many_dft_r2c(1,&rSize_[0],rSizeLocal_[1] ,rData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  //fPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,FFTW_FORWARD,FFTW_ESTIMATE);
  //fPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_FORWARD,FFTW_ESTIMATE);
  //fPlan_k_real_ =  fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE);

  //bPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],kSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,FFTW_BACKWARD,FFTW_ESTIMATE);
  //bPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_BACKWARD,FFTW_ESTIMATE);
  //bPlan_j_real_ =  fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE);
  bPlan_i_ = fftwf_plan_many_dft_c2r(1,&rSize_[0],rSizeLocal_[1] ,temp1_,NULL,1, r2cSize_,rData_,NULL,components_,rJump_[1]*components_,FFTW_ESTIMATE);

//  PrintPlans();

  //allocation of field

  long rfield_size = rfield->lattice().sitesLocalGross();
  long kfield_size = kfield->lattice().sitesLocalGross();

  if(mem_type_==FFT_IN_PLACE)
  {
    if(rfield_size>kfield_size*2)
    {
      rfield->alloc();
      kfield->data() = (Imag *)rfield->data();
    }
    else
    {
      kfield->alloc();
      rfield->data() = (float *)kfield->data();
    }
  }
  if(mem_type_ == FFT_OUT_OF_PLACE)
  {
    rfield->alloc(rfield_size, Field<float>::managed);
    kfield->alloc(kfield_size);
  }


  //Pointer to data

  rData_ = rfield->data();
  rData_ += rfield->lattice().siteFirst()*components_;//usefull point !!
  cData_ = (fftwf_complex*)kfield->data(); //to be sure that cData is instantiate !
  kData_ = (fftwf_complex*)kfield->data();
  kData_ += kfield->lattice().siteFirst()*components_;

  bool local_alignment_even = true;
#ifdef SINGLE
  if (get_pointer_alignment(rData_) < alignof(cufftComplex))
  {
    local_alignment_even = false;
  }
#else
  if (get_pointer_alignment(rData_) < alignof(cufftDoubleComplex))
  {
    local_alignment_even = false;
  }
#endif

  int local_alignment = local_alignment_even ? 1 : 0;
  int min_alignment = 0;
  int max_alignment = 0;
  MPI_Allreduce(&local_alignment, &min_alignment, 1, MPI_INT, MPI_MIN, parallel.lat_world_comm());
  MPI_Allreduce(&local_alignment, &max_alignment, 1, MPI_INT, MPI_MAX, parallel.lat_world_comm());

  if (min_alignment != max_alignment)
  {
    if(parallel.isRoot())
    {
      cerr<<"Latfield2d::PlanFFT::initialize : inconsistent data alignment across MPI processes for CUDA FFT"<<endl;
      cerr<<"Latfield2d : Abort Process Requested"<<endl;
    }
    parallel.abortForce();
  }

  alignment_even_ = (max_alignment == 1);



}

#endif

#ifndef SINGLE

template<class compType>
PlanFFT<compType>::PlanFFT(Field<compType>*  rfield,Field<compType>* kfield,const int mem_type, const bool managed)
{
  managed_ = managed;
	status_ = false;
	initialize(rfield,kfield,mem_type);
}

template<class compType>
void PlanFFT<compType>::initialize(Field<compType>*  rfield,Field<compType>*  kfield,const int mem_type )
{
  type_ = C2C;
  mem_type_=mem_type;

  //general variable

  if(rfield->components() != kfield->components())
  {
    cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for fields with same number of components"<<endl;
    cerr<<"Latfield2d::PlanFFT::initialize : coordinate and fourier space fields have not the same number of components"<<endl;
    cerr<<"Latfield2d : Abort Process Requested"<<endl;

  }
  else components_ = rfield->components();

  for(int i = 0; i<3; i++)
  {
    rSize_[i]=rfield->lattice().size(i);
    kSize_[i]=kfield->lattice().size(i);
    rSizeLocal_[i]=rfield->lattice().sizeLocal(i);
    kSizeLocal_[i]=kfield->lattice().sizeLocal(i);
    rJump_[i]=rfield->lattice().jump(i);
    kJump_[i]=kfield->lattice().jump(i);
  }
  r2cSize_ = 0;
  r2cSizeLocal_as_ = 0;
  r2cSizeLocal_ = 0;
  rHalo_ = rfield->lattice().halo();
  kHalo_ = kfield->lattice().halo();

  /////from latfield2d_IO
  tempMemory.setTemp((long)(rSize_[0] + 2*rHalo_)  * (long)(rSizeLocal_[1] + 2*rHalo_) * (long)(rSizeLocal_[2] + 2*rHalo_));

  	temp_  = tempMemory.temp1();
  	temp1_ = tempMemory.temp2();

  	if(rfield->lattice().dim()!=3)
  	{
  		if(parallel.isRoot())
  		{
  			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
  			cerr<<"Latfield2d::PlanFFT::initialize : real lattice have not 3 dimensions"<<endl;
  			cerr<<"Latfield2d : Abort Process Requested"<<endl;

  		}
  		parallel.abortForce();
  	}


  	if(rSize_[0]!=rSize_[1] | rSize_[1]!=rSize_[2])
  	{
  		if(parallel.isRoot())
  		{
  			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
  			cerr<<"Latfield2d::PlanFFT::initialize : real lattice is not cubic"<<endl;
  			cerr<<"Latfield2d : Abort Process Requested"<<endl;

   		}
  		parallel.abortForce();
  	}

  	//initialization of fftw plan


  	//Forward plan
  	fPlan_i_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,cData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  	fPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE);
  	fPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,kData_,NULL,components_, rJump_[1]*components_,FFTW_FORWARD,FFTW_ESTIMATE);
  	//Backward plan

  	bPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,kData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  	bPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE);
  	bPlan_i_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,cData_,NULL,components_, rJump_[1]*components_,FFTW_BACKWARD,FFTW_ESTIMATE);

  	//allocation of field

  	long rfield_size = rfield->lattice().sitesLocalGross();
  	long kfield_size = kfield->lattice().sitesLocalGross();

  	if(mem_type_ == FFT_IN_PLACE)
  	{
  		if(rfield_size>=kfield_size)
  		{
  			rfield->alloc();
  			kfield->data() = (Imag*)rfield->data();
  		}
  		else
  		{
  			kfield->alloc();
  			rfield->data() = (Imag*)kfield->data();
  		}
  	}
  	if(mem_type_ == FFT_OUT_OF_PLACE)
  	{
  		rfield->alloc();
  		kfield->alloc();
  	}

  	//Pointer to data

  	rData_ = (double*)rfield->data(); //to be sure that rData is instantiate !
  	cData_ = (fftw_complex*)rfield->data();
  	cData_ += rfield->lattice().siteFirst()*components_;
  	kData_ = (fftw_complex*)kfield->data();
  	kData_ += kfield->lattice().siteFirst()*components_;

}

template<class compType>
PlanFFT<compType>::PlanFFT(Field<double>* rfield, Field<compType>*  kfield,const int mem_type, const bool managed )
{
  managed_ = managed;
  status_ = false;
  initialize(rfield,kfield,mem_type);
}



template<class compType>
void PlanFFT<compType>::initialize(Field<double>*  rfield,Field<compType>*   kfield,const int mem_type )
{
  type_ = R2C;
  mem_type_=mem_type;

  //general variable
  if(rfield->components() != kfield->components())
  {
    cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for fields with same number of components"<<endl;
    cerr<<"Latfield2d::PlanFFT::initialize : coordinate and fourier space fields have not the same number of components"<<endl;
    cerr<<"Latfield2d : Abort Process Requested"<<endl;

  }
  else components_ = rfield->components();

  for(int i = 0; i<3; i++)
  {
    rSize_[i]=rfield->lattice().size(i);
    kSize_[i]=kfield->lattice().size(i);
    rSizeLocal_[i]=rfield->lattice().sizeLocal(i);
    kSizeLocal_[i]=kfield->lattice().sizeLocal(i);
    rJump_[i]=rfield->lattice().jump(i);
    kJump_[i]=kfield->lattice().jump(i);
  }
  r2cSize_ = rfield->lattice().size(0)/2 + 1;
  r2cSizeLocal_as_ = rfield->lattice().sizeLocal(0)/(2*parallel.grid_size()[1]);
  if(parallel.last_proc()[1]) r2cSizeLocal_ = r2cSizeLocal_as_ + 1;
  else r2cSizeLocal_ = r2cSizeLocal_as_;
  rHalo_ = rfield->lattice().halo();
  kHalo_ = kfield->lattice().halo();

  preallocate();

  temp_  = tempMemory.temp1();
  temp1_ = tempMemory.temp2();

  if(rfield->lattice().dim()!=3)
  {
    if(parallel.isRoot())
    {
      cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
      cerr<<"Latfield2d::PlanFFT::initialize : real lattice have not 3 dimensions"<<endl;
      cerr<<"Latfield2d : Abort Process Requested"<<endl;

    }
    parallel.abortForce();
  }

  //COUT << rSize_[0] << "   "<< rSize_[1]<< "   "<< rSize_[2] <<endl;

  if(rSize_[0]!=rSize_[1] | rSize_[1]!=rSize_[2])
  {
    if(parallel.isRoot())
    {
      cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
      cerr<<"Latfield2d::PlanFFT::initialize : real lattice is not cubic"<<endl;
      cerr<<"Latfield2d : Abort Process Requested"<<endl;

    }
    parallel.abortForce();
  }

  //create the fftw_plan

  long rfield_size = rfield->lattice().sitesLocalGross();
  long kfield_size = kfield->lattice().sitesLocalGross();

  if(mem_type_==FFT_IN_PLACE)
  {
    if(rfield_size>kfield_size*2)
    {
      rfield->alloc();
      kfield->data() = (Imag *)rfield->data();
    }
    else
    {
      kfield->alloc();
      rfield->data() = (double *)kfield->data();
    }
  }
  if(mem_type_ == FFT_OUT_OF_PLACE)
  {
    rfield->alloc(rfield_size, Field<double>::managed);
    kfield->alloc(kfield_size);
  }


  //Pointer to data

  rData_ = rfield->data();
  rData_ += rfield->lattice().siteFirst()*components_;//usefull point !!
  cData_ = (fftw_complex*)kfield->data(); //to be sure that cData is instantiate !
  kData_ = (fftw_complex*)kfield->data();
  kData_ += kfield->lattice().siteFirst()*components_;

  bool local_alignment_even = true;
#ifdef SINGLE
  if (get_pointer_alignment(rData_) < alignof(cufftComplex))
  {
    local_alignment_even = false;
  }
#else
  if (get_pointer_alignment(rData_) < alignof(cufftDoubleComplex))
  {
    local_alignment_even = false;
  }
#endif

  int local_alignment = local_alignment_even ? 1 : 0;
  int min_alignment = 0;
  int max_alignment = 0;
  MPI_Allreduce(&local_alignment, &min_alignment, 1, MPI_INT, MPI_MIN, parallel.lat_world_comm());
  MPI_Allreduce(&local_alignment, &max_alignment, 1, MPI_INT, MPI_MAX, parallel.lat_world_comm());

  if (min_alignment != max_alignment)
  {
    if(parallel.isRoot())
    {
      cerr<<"Latfield2d::PlanFFT::initialize : inconsistent data alignment across MPI processes for CUDA FFT"<<endl;
      cerr<<"Latfield2d : Abort Process Requested"<<endl;
    }
    parallel.abortForce();
  }

  alignment_even_ = (max_alignment == 1);

  cufftResult cufft_status;
  int inembed[1] = {rSizeLocal_[0]};
  int onembed[1] = {(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2] + 1};

  cufft_status = cufftPlanMany(&cufPlan_i_, 1, &rSize_[0], inembed, components_, rJump_[1]*components_, onembed, rSizeLocal_[1]*rSizeLocal_[2], 1, CUFFT_D2Z, rSizeLocal_[1]);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for forward plan i" << endl;
    parallel.abortForce();
  }

  inembed[0] = 1;
  onembed[0] = 1;

  cufft_status = cufftPlanMany(&cufPlan_j_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, 1, onembed, rSizeLocal_[2]*r2cSizeLocal_, 1, CUFFT_Z2Z, rSizeLocal_[2]*r2cSizeLocal_);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for forward plan j" << endl;
    parallel.abortForce();
  }

  cufft_status = cufftPlanMany(&cufPlan_k_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, 1, onembed, rSizeLocal_[2]*r2cSizeLocal_as_, 1, CUFFT_Z2Z, r2cSizeLocal_as_);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for forward plan k" << endl;
    parallel.abortForce();
  }

  inembed[0] = r2cSizeLocal_;

  cufft_status = cufftPlanMany(&cufPlan_k_real_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, r2cSizeLocal_, onembed, rSizeLocal_[2], 1, CUFFT_Z2Z, rSizeLocal_[2]);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for forward plan k_real" << endl;
    parallel.abortForce();
  }

  inembed[0] = 1;

  cufft_status = cufftPlanMany(&cubPlan_k_, 1, &rSize_[0], inembed, kSizeLocal_[2]*r2cSizeLocal_, 1, onembed, kSizeLocal_[2]*r2cSizeLocal_, 1, CUFFT_Z2Z, kSizeLocal_[2]*r2cSizeLocal_);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for backward plan k" << endl;
    parallel.abortForce();
  }

  cufft_status = cufftPlanMany(&cubPlan_j_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, 1, onembed, rSizeLocal_[2]*r2cSizeLocal_as_, 1, CUFFT_Z2Z, r2cSizeLocal_as_);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for backward plan j" << endl;
    parallel.abortForce();
  }

  inembed[0] = r2cSizeLocal_;

  cufft_status = cufftPlanMany(&cubPlan_j_real_, 1, &rSize_[0], inembed, rSizeLocal_[2]*r2cSizeLocal_, r2cSizeLocal_, onembed, rSizeLocal_[2], 1, CUFFT_Z2Z, rSizeLocal_[2]);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for backward plan j_real" << endl;
    parallel.abortForce();
  }

  inembed[0] = r2cSize_;
  onembed[0] = rJump_[1];

  cufft_status = cufftPlanMany(&cubPlan_i_, 1, &rSize_[0], inembed, 1, r2cSize_, onembed, components_, rJump_[1]*components_, CUFFT_Z2D, rSizeLocal_[1]);

  if(cufft_status != CUFFT_SUCCESS)
  {
    cerr << "Latfield2d::PlanFFT::initialize : cufftPlanMany failed for backward plan i" << endl;
    parallel.abortForce();
  }

  fPlan_i_ = fftw_plan_many_dft_r2c(1,&rSize_[0],rSizeLocal_[1] ,rData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  //fPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,FFTW_FORWARD,FFTW_ESTIMATE);
  //fPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_FORWARD,FFTW_ESTIMATE);
  //fPlan_k_real_ =  fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE);

  //bPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],kSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,FFTW_BACKWARD,FFTW_ESTIMATE);
  //bPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_BACKWARD,FFTW_ESTIMATE);
  //bPlan_j_real_ =  fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE);
  bPlan_i_ = fftw_plan_many_dft_c2r(1,&rSize_[0],rSizeLocal_[1] ,temp1_,NULL,1, r2cSize_,rData_,NULL,components_,rJump_[1]*components_,FFTW_ESTIMATE);

//  PrintPlans();
}

#endif

template<class compType>
void PlanFFT<compType>::preallocate()
{
  preallocate(0);
}

template<class compType>
long PlanFFT<compType>::requiredTemporaryCapacity() const
{
  long halo_pad = 2 * ((rHalo_ > kHalo_) ? rHalo_ : kHalo_);
  long required_real_dim0 = rSize_[0] + halo_pad;
  if (type_ == R2C) required_real_dim0 = (required_real_dim0 + 1) / 2;
  long required_real = required_real_dim0
                     * (long)(rSizeLocal_[1] + halo_pad)
                     * (long)(rSizeLocal_[2] + halo_pad);
  long required_fourier = (long)(((r2cSize_ > 0) ? r2cSize_ : rSize_[0]) + halo_pad)
                        * (long)(rSizeLocal_[1] + halo_pad)
                        * (long)(rSizeLocal_[2] + halo_pad);
  long required_capacity = required_real;
  if (required_fourier > required_capacity) required_capacity = required_fourier;

  if (type_ == R2C)
  {
    // The last dim1 rank keeps one extra Nyquist rod after the first gather.
    long required_redistributed = (long)parallel.grid_size()[1]
                                * (long)r2cSizeLocal_
                                * (long)rSizeLocal_[1]
                                * (long)rSizeLocal_[2];
    if (required_redistributed > required_capacity) required_capacity = required_redistributed;
  }

  return required_capacity;
}

template<class compType>
void PlanFFT<compType>::preallocate(size_t min_shared_device_workspace_bytes)
{
  tempMemory.setTemp(requiredTemporaryCapacity());
  tempMemory.reserveDeviceWorkspaceBytes(min_shared_device_workspace_bytes);
}

template<class compType>
void PlanFFT<compType>::execute(int fft_type)
{
  bool alignment_even = true;
  cudaStream_t fft_stream;

  long required_capacity = requiredTemporaryCapacity();

  if (tempMemory.allocated() < required_capacity)
  {
    // tempMemory.setTemp(required_capacity);
    COUT << "NB: required reallocation of memory! Should only be allocated once. Check initial allocation\n";

    if (tempMemory.allocated() < required_capacity)
    {
      cerr << "Latfield2d::PlanFFT::execute : insufficient temporary FFT buffer on proc "
           << parallel.rank() << " (required=" << required_capacity
           << ", allocated=" << tempMemory.allocated()
           << ", halo_real=" << rHalo_ << ", halo_fourier=" << kHalo_
           << ", local_size=(" << rSizeLocal_[0] << "," << rSizeLocal_[1] << "," << rSizeLocal_[2] << ")"
           << ", grid=" << parallel.grid_size()[0] << "x" << parallel.grid_size()[1] << ")"
           << endl;
      parallel.abortForce();
    }
  }

  temp_  = tempMemory.temp1();
  temp1_ = tempMemory.temp2();
  resolveCudaAwareMPI();
  const bool cuda_aware_mpi_active = cuda_aware_mpi_enabled_;

  auto success = cudaStreamCreateWithFlags(&fft_stream, cudaStreamNonBlocking);
  if (success != cudaSuccess)
  {
    cerr << "Latfield2d::PlanFFT::execute : cudaStreamCreateWithFlags failed" << endl;
    parallel.abortForce();
  }

  auto sync_fft_stream = [&](const char* context)
  {
    auto stream_status = cudaStreamSynchronize(fft_stream);
    if (stream_status != cudaSuccess)
    {
      cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed";
      if (context != nullptr) cerr << " during " << context;
      cerr << endl;
      parallel.abortForce();
    }
  };

  auto copy_device_buffer = [&](void* dst, const void* src, long bytes, const char* label)
  {
    if (label != nullptr) nvtxRangePushA(label);
    auto copy_status = cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToDevice, fft_stream);
    if (copy_status != cudaSuccess)
    {
      cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << copy_status << " (" << cudaGetErrorString((cudaError_t)copy_status) << ")" << endl;
      parallel.abortForce();
    }
    if (label != nullptr) nvtxRangePop();
  };

  auto complex_offset_ptr = [](void* ptr, long complex_offset) -> void*
  {
    return reinterpret_cast<void*>(reinterpret_cast<Real*>(ptr) + 2 * complex_offset);
  };

  auto complex_offset_ptr_const = [](const void* ptr, long complex_offset) -> const void*
  {
    return reinterpret_cast<const void*>(reinterpret_cast<const Real*>(ptr) + 2 * complex_offset);
  };

  auto mpi_alltoall_gather_cuda = [&](const char* label,
                                      bool barrier_world,
                                      const void* sendbuf,
                                      void* recvbuf,
                                      int alltoall_count,
                                      long gather_send_complex_offset,
                                      long gather_recv_complex_offset,
                                      int gather_count,
                                      MPI_Comm comm)
  {
    if (label != nullptr) nvtxRangePushA(label);
    sync_fft_stream(label);
    if (barrier_world) MPI_Barrier(parallel.lat_world_comm());
    MPI_Alltoall(const_cast<void*>(sendbuf), alltoall_count, MPI_DATA_PREC, recvbuf, alltoall_count, MPI_DATA_PREC, comm);
    MPI_Gather(const_cast<void*>(complex_offset_ptr_const(sendbuf, gather_send_complex_offset)), gather_count, MPI_DATA_PREC,
               complex_offset_ptr(recvbuf, gather_recv_complex_offset), gather_count, MPI_DATA_PREC,
               parallel.grid_size()[1]-1, comm);
    MPI_Barrier(comm);
    if (label != nullptr) nvtxRangePop();
  };

  auto mpi_alltoall_cuda = [&](const char* label,
                               bool barrier_world,
                               const void* sendbuf,
                               void* recvbuf,
                               int alltoall_count,
                               MPI_Comm comm)
  {
    if (label != nullptr) nvtxRangePushA(label);
    sync_fft_stream(label);
    if (barrier_world) MPI_Barrier(parallel.lat_world_comm());
    MPI_Alltoall(const_cast<void*>(sendbuf), alltoall_count, MPI_DATA_PREC, recvbuf, alltoall_count, MPI_DATA_PREC, comm);
    MPI_Barrier(comm);
    if (label != nullptr) nvtxRangePop();
  };

  auto mpi_alltoall_scatter_cuda = [&](const char* label,
                                       const void* sendbuf,
                                       void* recvbuf,
                                       int alltoall_count,
                                       long scatter_send_complex_offset,
                                       long scatter_recv_complex_offset,
                                       int scatter_count,
                                       MPI_Comm comm)
  {
    if (label != nullptr) nvtxRangePushA(label);
    sync_fft_stream(label);
    MPI_Alltoall(const_cast<void*>(sendbuf), alltoall_count, MPI_DATA_PREC, recvbuf, alltoall_count, MPI_DATA_PREC, comm);
    MPI_Scatter(const_cast<void*>(complex_offset_ptr_const(sendbuf, scatter_send_complex_offset)), scatter_count, MPI_DATA_PREC,
                complex_offset_ptr(recvbuf, scatter_recv_complex_offset), scatter_count, MPI_DATA_PREC,
                parallel.grid_size()[1]-1, comm);
    MPI_Barrier(comm);
    if (label != nullptr) nvtxRangePop();
  };

  //#ifdef SINGLE
  if(type_ == R2C)
  {

    if(fft_type == FFT_FORWARD)
    {
      //int i,j,k;
      //int comp;
      //int comm_rank;

#ifdef SINGLE
      float *p_in;
      fftwf_complex *p_out;
      cufftComplex *cutemp_ = tempMemory.temp3();
      cufftComplex *cutemp1_ = tempMemory.temp4();
      cufftComplex *cucomm_send_ = tempMemory.temp5();
      cufftComplex *p_cufft;
      long memsize = tempMemory.allocated()*sizeof(cufftComplex);
#else
      double *p_in;
      fftw_complex *p_out;
      cufftDoubleComplex *cutemp_ = tempMemory.temp3();
      cufftDoubleComplex *cutemp1_ = tempMemory.temp4();
      cufftDoubleComplex *cucomm_send_ = tempMemory.temp5();
      cufftDoubleComplex *p_cufft;
      long memsize = tempMemory.allocated()*sizeof(cufftDoubleComplex);
#endif
      alignment_even = alignment_even_;

      auto cufft_status = cufftSetStream(cufPlan_i_, fft_stream);
      if (cufft_status != CUFFT_SUCCESS)
      {
        cerr << "Latfield2d::PlanFFT::execute : cufftSetStream failed for plan i" << endl;
        parallel.abortForce();
      }

      cufft_status = cufftSetStream(cufPlan_j_, fft_stream);
      if (cufft_status != CUFFT_SUCCESS)
      {
        cerr << "Latfield2d::PlanFFT::execute : cufftSetStream failed for plan j" << endl;
        parallel.abortForce();
      }

      cufft_status = cufftSetStream(cufPlan_k_, fft_stream);
      if (cufft_status != CUFFT_SUCCESS)
      {
        cerr << "Latfield2d::PlanFFT::execute : cufftSetStream failed for plan k" << endl;
        parallel.abortForce();
      }

      cufft_status = cufftSetStream(cufPlan_k_real_, fft_stream);
      if (cufft_status != CUFFT_SUCCESS)
      {
        cerr << "Latfield2d::PlanFFT::execute : cufftSetStream failed for plan k_real" << endl;
        parallel.abortForce();
      }

      if (alignment_even || components_ > 1)
      {
        for(int l = 0; l < rSizeLocal_[2]; l++)
        {
          p_in = &rData_[rJump_[2]*l*components_ + (alignment_even ? 0 : 1)];
          p_cufft = &cutemp_[l*rSizeLocal_[1]];

  #ifdef SINGLE
          cufft_status = cufftExecR2C(cufPlan_i_, (cufftReal*)p_in, p_cufft);
  #else
          cufft_status = cufftExecD2Z(cufPlan_i_, (cufftDoubleReal*)p_in, p_cufft);
  #endif
          if (cufft_status != CUFFT_SUCCESS)
          {
  #ifdef SINGLE
            cerr << "Latfield2d::PlanFFT::execute : cufftExecR2C failed on component " << (alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
            cerr << "Latfield2d::PlanFFT::execute : cufftExecD2Z failed on component " << (alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
            parallel.abortForce();
          }
        }
      }

      for(int comp=0; comp<components_; comp += 2)
      {
        //execute first dimension fft, prepar rData in temp_ to be send via AlltoAll + gather

        nvtxRangePushA("PlanFFT::execute : R2C transforms");

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("PlanFFT::execute : interleaved R2C transform");
#pragma omp parallel for private(p_in, p_out)
          for(int l = 0; l < rSizeLocal_[2]; l++)
          {
            p_in = &rData_[rJump_[2]*l*components_ + comp + (alignment_even ? 1 : 0)];
            p_out = &temp_[l*rSizeLocal_[1]];
#ifdef SINGLE
            fftwf_execute_dft_r2c(fPlan_i_,p_in,p_out);
#else
            fftw_execute_dft_r2c(fPlan_i_,p_in,p_out);
#endif
          }
          nvtxRangePop();

          nvtxRangePushA("PlanFFT::execute : interleaved MPI AlltoAll and Gather");
          MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
          MPI_Gather(&temp_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
          MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
          nvtxRangePop();

          success = cudaMemcpyAsync(cutemp1_, temp1_, memsize, cudaMemcpyDefault, fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
            parallel.abortForce();
          }
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            copy_device_buffer(cucomm_send_, cutemp_, memsize, "Latfield2d::PlanFFT::execute : CUDA-aware MPI send staging");
          }
          else
          {
            success = cudaMemcpyAsync(temp_, cutemp_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }
        }

        success = cudaStreamSynchronize(fft_stream);
        if (success != cudaSuccess)
        {
          cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
          parallel.abortForce();
        }        
        nvtxRangePop();

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("PlanFFT::execute : interleaved transpose");
          if(parallel.last_proc()[1])
          {
            for(int i=0;i<parallel.grid_size()[1];i++)
            {
              //transpose_0_2_last_proc(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_);
              transpose_0_2_last_proc(&cutemp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&cutemp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_, fft_stream);
            }
            //implement_local_0_last_proc(&temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]],temp_,rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_,parallel.grid_size()[1]);
            implement_local_0_last_proc(&cutemp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]],cutemp_,rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_,parallel.grid_size()[1], fft_stream);
          }
          else
          {
            for(int i=0;i<parallel.grid_size()[1];i++)
            {
              //transpose_0_2(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_);
              transpose_0_2(&cutemp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&cutemp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_, fft_stream);
            }
          }
          nvtxRangePop();

          nvtxRangePushA("PlanFFT::execute : interleaved C2C transform");
  #ifdef SINGLE
          cufft_status = cufftExecC2C(cufPlan_j_, cutemp_, cutemp_, CUFFT_FORWARD);
  #else
          cufft_status = cufftExecZ2Z(cufPlan_j_, cutemp_, cutemp_, CUFFT_FORWARD);
  #endif
          if (cufft_status != CUFFT_SUCCESS)
          {
  #ifdef SINGLE
            cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
            cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
            parallel.abortForce();
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            mpi_alltoall_gather_cuda("PlanFFT::execute : CUDA-aware MPI AlltoAll and Gather",
                                     false,
                                     cucomm_send_,
                                     cutemp1_,
                                     2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_,
                                     (long)(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2],
                                     (long)(rSize_[0]/2)*rSizeLocal_[1]*rSizeLocal_[2],
                                     2*rSizeLocal_[1]*rSizeLocal_[2],
                                     parallel.dim1_comm()[parallel.grid_rank()[0]]);
          }
          else
          {
            nvtxRangePushA("PlanFFT::execute : MPI AlltoAll and Gather");
            MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Gather(&temp_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
            nvtxRangePop();

            success = cudaMemcpyAsync(cutemp1_, temp1_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("PlanFFT::execute : extra cudaMemcpyAsync due to interleaving");
          success = cudaMemcpyAsync(temp_, cutemp_, memsize, cudaMemcpyDefault, fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
            parallel.abortForce();
          }

          success = cudaStreamSynchronize(fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
            parallel.abortForce();
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || alignment_even)
        {
          nvtxRangePushA("PlanFFT::execute : Transpose");
          if(parallel.last_proc()[1])
          {
            for(int i=0;i<parallel.grid_size()[1];i++)
            {
              transpose_0_2_last_proc(&cutemp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&cutemp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_, fft_stream);
            }
            implement_local_0_last_proc(&cutemp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]],cutemp_,rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_,parallel.grid_size()[1], fft_stream);
          }
          else
          {
            for(int i=0;i<parallel.grid_size()[1];i++)
            {
              transpose_0_2(&cutemp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&cutemp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_, fft_stream);
            }
          }
          nvtxRangePop();

          nvtxRangePushA("PlanFFT::execute : C2C transform");
  #ifdef SINGLE
          cufft_status = cufftExecC2C(cufPlan_j_, cutemp_, cutemp_, CUFFT_FORWARD);
  #else
          cufft_status = cufftExecZ2Z(cufPlan_j_, cutemp_, cutemp_, CUFFT_FORWARD);
  #endif
          if (cufft_status != CUFFT_SUCCESS)
          {
  #ifdef SINGLE
            cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
            cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
            parallel.abortForce();
          }
        }
        
        //fftwf_execute(fPlan_j_);

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("PlanFFT::execute : interleaved MPI AlltoAll");
          MPI_Barrier(parallel.lat_world_comm());
          MPI_Alltoall(temp_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, temp1_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);
          MPI_Barrier(parallel.dim0_comm()[parallel.grid_rank()[1]]);
          nvtxRangePop();

          success = cudaMemcpyAsync(cutemp1_, temp1_, memsize, cudaMemcpyDefault, fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
            parallel.abortForce();
          }
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            copy_device_buffer(cucomm_send_, cutemp_, memsize, "PlanFFT::execute : CUDA-aware MPI send staging");
          }
          else
          {
            success = cudaMemcpyAsync(temp_, cutemp_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }
          nvtxRangePop();
        }

        success = cudaStreamSynchronize(fft_stream);
        if (success != cudaSuccess)
        {
          cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
          parallel.abortForce();
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("PlanFFT::execute : interleaved transpose");
          for(int i=0;i<parallel.grid_size()[0];i++)
          {
            //transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_], r2cSizeLocal_,rSizeLocal_[2],rSizeLocal_[2]);
            transpose_1_2(&cutemp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_],&cutemp_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_], r2cSizeLocal_,rSizeLocal_[2],rSizeLocal_[2], fft_stream);
          }
          nvtxRangePop();

          nvtxRangePushA("PlanFFT::execute : interleaved C2C transform");
          for(int l=0;l<rSizeLocal_[2];l++)
          {
            //fftwf_execute_dft(fPlan_k_,&temp_[l*r2cSizeLocal_],&temp1_[l*r2cSizeLocal_as_]);
  #ifdef SINGLE
            cufft_status = cufftExecC2C(cufPlan_k_, &cutemp_[l*r2cSizeLocal_], &cutemp1_[l*r2cSizeLocal_as_], CUFFT_FORWARD);
  #else
            cufft_status = cufftExecZ2Z(cufPlan_k_, &cutemp_[l*r2cSizeLocal_], &cutemp1_[l*r2cSizeLocal_as_], CUFFT_FORWARD);
  #endif
            if (cufft_status != CUFFT_SUCCESS)
            {
  #ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
              parallel.abortForce();
            }
          }

          if(parallel.last_proc()[1]) //fftwf_execute(fPlan_k_real_);
          {
  #ifdef SINGLE
            cufft_status = cufftExecC2C(cufPlan_k_real_, &cutemp_[r2cSizeLocal_as_], &cutemp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]], CUFFT_FORWARD);
  #else
            cufft_status = cufftExecZ2Z(cufPlan_k_real_, &cutemp_[r2cSizeLocal_as_], &cutemp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]], CUFFT_FORWARD);
  #endif
            if (cufft_status != CUFFT_SUCCESS)
            {
  #ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
              parallel.abortForce();
            }
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            mpi_alltoall_cuda("PlanFFT::execute : CUDA-aware MPI AlltoAll",
                              true,
                              cucomm_send_,
                              cutemp_,
                              (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_),
                              parallel.dim0_comm()[parallel.grid_rank()[1]]);
          }
          else
          {
            nvtxRangePushA("PlanFFT::execute : MPI AlltoAll");
            MPI_Barrier(parallel.lat_world_comm());
            MPI_Alltoall(temp_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, temp1_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);
            MPI_Barrier(parallel.dim0_comm()[parallel.grid_rank()[1]]);
            nvtxRangePop();

            success = cudaMemcpyAsync(cutemp_, temp1_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("PlanFFT::execute : extra cudaMemcpyAsync due to interleaving");
          success = cudaMemcpyAsync(temp_, cutemp1_, memsize, cudaMemcpyDefault, fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
            parallel.abortForce();
          }

          success = cudaStreamSynchronize(fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
            parallel.abortForce();
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || alignment_even)
        {
          nvtxRangePushA("PlanFFT::execute : Transpose");
          for(int i=0;i<parallel.grid_size()[0];i++)
          {
            transpose_1_2(&cutemp_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_],&cutemp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_], r2cSizeLocal_,rSizeLocal_[2],rSizeLocal_[2], fft_stream);
          }
          nvtxRangePop();

          nvtxRangePushA("PlanFFT::execute : C2C transform");
          for(int l=0;l<rSizeLocal_[2];l++)
          {
            //fftwf_execute_dft(fPlan_k_,&temp_[l*r2cSizeLocal_],&temp1_[l*r2cSizeLocal_as_]);
  #ifdef SINGLE
            cufft_status = cufftExecC2C(cufPlan_k_, &cutemp1_[l*r2cSizeLocal_], &cutemp_[l*r2cSizeLocal_as_], CUFFT_FORWARD);
  #else
            cufft_status = cufftExecZ2Z(cufPlan_k_, &cutemp1_[l*r2cSizeLocal_], &cutemp_[l*r2cSizeLocal_as_], CUFFT_FORWARD);
  #endif
            if (cufft_status != CUFFT_SUCCESS)
            {
  #ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+1 << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+1 << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
              parallel.abortForce();
            }
          }

          if(parallel.last_proc()[1]) //fftwf_execute(fPlan_k_real_);
          {
  #ifdef SINGLE
            cufft_status = cufftExecC2C(cufPlan_k_real_, &cutemp1_[r2cSizeLocal_as_], &cutemp_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]], CUFFT_FORWARD);
  #else
            cufft_status = cufftExecZ2Z(cufPlan_k_real_, &cutemp1_[r2cSizeLocal_as_], &cutemp_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]], CUFFT_FORWARD);
  #endif
            if (cufft_status != CUFFT_SUCCESS)
            {
  #ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
              parallel.abortForce();
            }
          }
          nvtxRangePop();
        }
        
        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("PlanFFT::execute : interleaved MPI AlltoAll and Scatter");
          MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
          MPI_Scatter(&temp_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp1_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
          MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
          nvtxRangePop();

          nvtxRangePushA("PlanFFT::execute : interleaved transpose back");
          transpose_back_0_3(temp1_, kData_,r2cSize_,r2cSizeLocal_as_,rSizeLocal_[2],rSizeLocal_[1],parallel.grid_size()[1],kHalo_,components_,comp+(alignment_even ? 1 : 0));
          implement_0(&temp1_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]], kData_,r2cSize_,rSizeLocal_[2],rSizeLocal_[1],kHalo_,components_,comp+(alignment_even ? 1 : 0));
          nvtxRangePop();
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            copy_device_buffer(cucomm_send_, cutemp_, memsize, "PlanFFT::execute : CUDA-aware MPI send staging");
          }
          else
          {
            success = cudaMemcpyAsync(temp1_, cutemp_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }
          sync_fft_stream(comp + 1 < components_ || alignment_even ? "forward active scatter staging" : nullptr);
        }

        if (comp + (alignment_even ? 2 : 3) < components_) // offload next component
        {
          for(int l = 0;l< rSizeLocal_[2] ;l++)
          {
            p_in = &rData_[rJump_[2]*l*components_ + comp + (alignment_even ? 2 : 3)];
            p_cufft = &cutemp_[l*rSizeLocal_[1]];
            
#ifdef SINGLE
            cufft_status = cufftExecR2C(cufPlan_i_, (cufftReal*)p_in, p_cufft);
#else
            cufft_status = cufftExecD2Z(cufPlan_i_, (cufftDoubleReal*)p_in, p_cufft);
#endif
            if (cufft_status != CUFFT_SUCCESS)
            {
#ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecR2C failed on component " << comp+(alignment_even ? 3 : 4) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
#else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecD2Z failed on component " << comp+(alignment_even ? 3 : 4) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
#endif
              parallel.abortForce();
            }
          }
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            mpi_alltoall_scatter_cuda("PlanFFT::execute : CUDA-aware MPI AlltoAll and Scatter",
                                      cucomm_send_,
                                      temp_,
                                      2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_,
                                      (long)r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0],
                                      (long)(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2],
                                      2*rSizeLocal_[1]*rSizeLocal_[2],
                                      parallel.dim1_comm()[parallel.grid_rank()[0]]);
          }
          else
          {
            nvtxRangePushA("PlanFFT::execute : MPI AlltoAll and Scatter");
            MPI_Alltoall(temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Scatter(&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
            nvtxRangePop();
          }

          nvtxRangePushA("PlanFFT::execute : Transpose back");
          transpose_back_0_3(temp_, kData_,r2cSize_,r2cSizeLocal_as_,rSizeLocal_[2],rSizeLocal_[1],parallel.grid_size()[1],kHalo_,components_,comp+(alignment_even ? 0 : 1));
          implement_0(&temp_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]], kData_,r2cSize_,rSizeLocal_[2],rSizeLocal_[1],kHalo_,components_,comp+(alignment_even ? 0 : 1));
          nvtxRangePop();
        }
      }
    }
    if(fft_type == FFT_BACKWARD)
    {
      //int i,j,k,comp;
      //int comm_rank;

#ifdef SINGLE
      float *p_out;
      fftwf_complex *p_in;
      cufftComplex *cutemp_ = tempMemory.temp3();
      cufftComplex *cutemp1_ = tempMemory.temp4();
      cufftComplex *cucomm_send_ = tempMemory.temp5();
      cufftComplex *p_cufft;
      long memsize = tempMemory.allocated()*sizeof(cufftComplex);
#else
      double *p_out;
      fftw_complex *p_in;
      cufftDoubleComplex *cutemp_ = tempMemory.temp3();
      cufftDoubleComplex *cutemp1_ = tempMemory.temp4();
      cufftDoubleComplex *cucomm_send_ = tempMemory.temp5();
      cufftDoubleComplex *p_cufft;
      long memsize = tempMemory.allocated()*sizeof(cufftDoubleComplex);
#endif
      alignment_even = alignment_even_;

      auto cufft_status = cufftSetStream(cubPlan_k_, fft_stream);
      if (cufft_status != CUFFT_SUCCESS)
      {
        cerr << "Latfield2d::PlanFFT::execute : cufftSetStream failed for plan k" << endl;
        parallel.abortForce();
      }

      cufft_status = cufftSetStream(cubPlan_j_, fft_stream);
      if (cufft_status != CUFFT_SUCCESS)
      {
        cerr << "Latfield2d::PlanFFT::execute : cufftSetStream failed for plan j" << endl;
        parallel.abortForce();
      }

      cufft_status = cufftSetStream(cubPlan_j_real_, fft_stream);
      if (cufft_status != CUFFT_SUCCESS)
      {
        cerr << "Latfield2d::PlanFFT::execute : cufftSetStream failed for plan j_real" << endl;
        parallel.abortForce();
      }

      cufft_status = cufftSetStream(cubPlan_i_, fft_stream);
      if (cufft_status != CUFFT_SUCCESS)
      {
        cerr << "Latfield2d::PlanFFT::execute : cufftSetStream failed for plan i_real" << endl;
        parallel.abortForce();
      }

      for(int comp=0; comp<components_; comp += 2)
      {
        if (comp + 1 < components_ || alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : MPI AlltoAll and Gather");
          b_arrange_data_0(kData_, temp_,kSizeLocal_[0],kSizeLocal_[1] ,kSizeLocal_[2], kHalo_, components_, comp + (alignment_even ? 0 : 1));
          MPI_Barrier(parallel.lat_world_comm());
          if (cuda_aware_mpi_active)
          {
            MPI_Alltoall(temp_,2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, cutemp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Gather(&temp_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC,
                       complex_offset_ptr(cutemp1_, (long)(rSize_[0]/2)*rSizeLocal_[1]*rSizeLocal_[2]), 2*rSizeLocal_[1]*rSizeLocal_[2],
                       MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
          }
          else
          {
            MPI_Alltoall(temp_,2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Gather(&temp_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
          }
          nvtxRangePop();
        }

        if (comp > 0) // deferred synchronisation
        {
          success = cudaStreamSynchronize(fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
            parallel.abortForce();
          }
        }

        if (comp + 1 < components_ || alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : Transpose");
          if (!cuda_aware_mpi_active)
          {
            success = cudaMemcpyAsync(cutemp1_, temp1_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }

          if(parallel.last_proc()[1])
          {
            for(int i=0;i<parallel.grid_size()[1];i++)
            {
              transpose_0_2_last_proc(&cutemp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&cutemp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_, fft_stream);
            }
            implement_local_0_last_proc(&cutemp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]],cutemp_,rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_,parallel.grid_size()[1], fft_stream);
          }
          else
          {
            for(int i=0;i<parallel.grid_size()[1];i++)
            {
              transpose_0_2(&cutemp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&cutemp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_, fft_stream);
            }
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved b_arrange_data_0");
          b_arrange_data_0(kData_, temp_,kSizeLocal_[0],kSizeLocal_[1] ,kSizeLocal_[2], kHalo_, components_, comp + (alignment_even ? 1 : 0));

          success = cudaStreamSynchronize(fft_stream); // extra stream synchronisation before interleaved MPI
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
            parallel.abortForce();
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : C2C transform");
  #ifdef SINGLE
          //fftwf_execute(bPlan_k_);
          cufft_status = cufftExecC2C(cubPlan_k_, cutemp_, cutemp_, CUFFT_INVERSE);
  #else
          cufft_status = cufftExecZ2Z(cubPlan_k_, cutemp_, cutemp_, CUFFT_INVERSE);
  #endif
          if (cufft_status != CUFFT_SUCCESS)
          {
  #ifdef SINGLE
            cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
            cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
            parallel.abortForce();
          }
          //fftwf_execute(fPlan_j_);
          nvtxRangePop();
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved MPI AlltoAll and Gather");
          MPI_Barrier(parallel.lat_world_comm());

          MPI_Alltoall(temp_,2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
          MPI_Gather(&temp_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
          MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);

          success = cudaMemcpyAsync(cutemp1_, temp1_, memsize, cudaMemcpyDefault, fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
            parallel.abortForce();
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            copy_device_buffer(cucomm_send_, cutemp_, memsize, "Latfield2d::PlanFFT::execute : CUDA-aware MPI send staging");
          }
          else
          {
            success = cudaMemcpyAsync(temp_, cutemp_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }
        }

        success = cudaStreamSynchronize(fft_stream);
        if (success != cudaSuccess)
        {
          cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
          parallel.abortForce();
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved transpose");
          if(parallel.last_proc()[1])
          {
            for(int i=0;i<parallel.grid_size()[1];i++)
            {
              //transpose_0_2_last_proc(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_);
              transpose_0_2_last_proc(&cutemp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&cutemp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_, fft_stream);
            }
            //implement_local_0_last_proc(&temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]],temp_,rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_,parallel.grid_size()[1]);
            implement_local_0_last_proc(&cutemp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]],cutemp_,rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_,parallel.grid_size()[1], fft_stream);
          }
          else
          {
            for(int i=0;i<parallel.grid_size()[1];i++)
            {
              //transpose_0_2(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_);
              transpose_0_2(&cutemp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&cutemp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_, fft_stream);
            }
          }
          nvtxRangePop();
          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved C2C transform");
  #ifdef SINGLE
          //fftwf_execute(bPlan_k_);
          cufft_status = cufftExecC2C(cubPlan_k_, cutemp_, cutemp_, CUFFT_INVERSE);
  #else
          cufft_status = cufftExecZ2Z(cubPlan_k_, cutemp_, cutemp_, CUFFT_INVERSE);
  #endif
          if (cufft_status != CUFFT_SUCCESS)
          {
  #ifdef SINGLE
            cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
            cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
            parallel.abortForce();
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            mpi_alltoall_cuda("Latfield2d::PlanFFT::execute : CUDA-aware MPI AlltoAll",
                              false,
                              cucomm_send_,
                              cutemp1_,
                              (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_),
                              parallel.dim0_comm()[parallel.grid_rank()[1]]);
          }
          else
          {
            nvtxRangePushA("Latfield2d::PlanFFT::execute : MPI AlltoAll");
            MPI_Alltoall(temp_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, temp1_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);
            MPI_Barrier(parallel.dim0_comm()[parallel.grid_rank()[1]]);
            nvtxRangePop();

            success = cudaMemcpyAsync(cutemp1_, temp1_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          success = cudaMemcpyAsync(temp_, cutemp_, memsize, cudaMemcpyDefault, fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
            parallel.abortForce();
          }

          success = cudaStreamSynchronize(fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
            parallel.abortForce();
          }
        }

        if (comp + 1 < components_ || alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : Transpose");
          for(int i=0;i<parallel.grid_size()[0];i++)
          {
            transpose_1_2(&cutemp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_],&cutemp_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_], r2cSizeLocal_,rSizeLocal_[2],rSizeLocal_[2], fft_stream);
          }
          nvtxRangePop();

          nvtxRangePushA("Latfield2d::PlanFFT::execute : C2C transform");
          for(int l=0;l<rSizeLocal_[2];l++)
          {
            //fftwf_execute_dft(bPlan_j_,&temp_[l*r2cSizeLocal_],&temp1_[l*r2cSizeLocal_as_]);
  #ifdef SINGLE
            cufft_status = cufftExecC2C(cubPlan_j_, &cutemp_[l*r2cSizeLocal_], &cutemp1_[l*r2cSizeLocal_as_], CUFFT_INVERSE);
  #else
            cufft_status = cufftExecZ2Z(cubPlan_j_, &cutemp_[l*r2cSizeLocal_], &cutemp1_[l*r2cSizeLocal_as_], CUFFT_INVERSE);
  #endif
            if (cufft_status != CUFFT_SUCCESS)
            {
  #ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
              parallel.abortForce();
            }
          }

          if(parallel.last_proc()[1])
          {
            //fftwf_execute(bPlan_j_real_);
  #ifdef SINGLE
            cufft_status = cufftExecC2C(cubPlan_j_real_, &cutemp_[r2cSizeLocal_as_], &cutemp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]], CUFFT_INVERSE);
  #else
            cufft_status = cufftExecZ2Z(cubPlan_j_real_, &cutemp_[r2cSizeLocal_as_], &cutemp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]], CUFFT_INVERSE);
  #endif
            if (cufft_status != CUFFT_SUCCESS)
            {
  #ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
              parallel.abortForce();
            }
          }
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved MPI AlltoAll");
          MPI_Alltoall(temp_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, temp1_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);
          MPI_Barrier(parallel.dim0_comm()[parallel.grid_rank()[1]]);
          nvtxRangePop();

          success = cudaMemcpyAsync(cutemp_, temp1_, memsize, cudaMemcpyDefault, fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
            parallel.abortForce();
          }
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            copy_device_buffer(cucomm_send_, cutemp1_, memsize, "Latfield2d::PlanFFT::execute : CUDA-aware MPI send staging");
          }
          else
          {
            success = cudaMemcpyAsync(temp_, cutemp1_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }
        }

        success = cudaStreamSynchronize(fft_stream);
        if (success != cudaSuccess)
        {
          cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
          parallel.abortForce();
        }
        nvtxRangePop();

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved transpose");
          for(int i=0;i<parallel.grid_size()[0];i++)
          {
            //transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_], r2cSizeLocal_,rSizeLocal_[2],rSizeLocal_[2]);
            transpose_1_2(&cutemp_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_],&cutemp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_], r2cSizeLocal_,rSizeLocal_[2],rSizeLocal_[2], fft_stream);
          }
          nvtxRangePop();

          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved C2C transform");
          for(int l=0;l<rSizeLocal_[2];l++)
          {
            //fftwf_execute_dft(bPlan_j_,&temp_[l*r2cSizeLocal_],&temp1_[l*r2cSizeLocal_as_]);
#ifdef SINGLE
            cufft_status = cufftExecC2C(cubPlan_j_, &cutemp1_[l*r2cSizeLocal_], &cutemp_[l*r2cSizeLocal_as_], CUFFT_INVERSE);
#else
            cufft_status = cufftExecZ2Z(cubPlan_j_, &cutemp1_[l*r2cSizeLocal_], &cutemp_[l*r2cSizeLocal_as_], CUFFT_INVERSE);
#endif
            if (cufft_status != CUFFT_SUCCESS)
            {
#ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+2 << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
#else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+2 << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
#endif
              parallel.abortForce();
            }
          }

          if(parallel.last_proc()[1])
          {
            //fftwf_execute(bPlan_j_real_);
#ifdef SINGLE
            cufft_status = cufftExecC2C(cubPlan_j_real_, &cutemp1_[r2cSizeLocal_as_], &cutemp_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]], CUFFT_INVERSE);
#else
            cufft_status = cufftExecZ2Z(cubPlan_j_real_, &cutemp1_[r2cSizeLocal_as_], &cutemp_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]], CUFFT_INVERSE);
#endif
            if (cufft_status != CUFFT_SUCCESS)
            {
#ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecC2C failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
#else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2Z failed on component " << comp+(alignment_even ? 2 : 1) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
#endif
              parallel.abortForce();
            }
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || alignment_even)
        {
          if (cuda_aware_mpi_active)
          {
            mpi_alltoall_scatter_cuda("Latfield2d::PlanFFT::execute : CUDA-aware MPI AlltoAll and Scatter",
                                      cucomm_send_,
                                      cutemp1_,
                                      2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_,
                                      (long)r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0],
                                      (long)(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2],
                                      2*rSizeLocal_[1]*rSizeLocal_[2],
                                      parallel.dim1_comm()[parallel.grid_rank()[0]]);
          }
          else
          {
            nvtxRangePushA("Latfield2d::PlanFFT::execute : MPI AlltoAll and Scatter");
            MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Scatter(&temp_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp1_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
            MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
            nvtxRangePop();

            success = cudaMemcpyAsync(cutemp1_, temp1_, memsize, cudaMemcpyDefault, fft_stream);
            if (success != cudaSuccess)
            {
              cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
              parallel.abortForce();
            }
          }
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          success = cudaMemcpyAsync(temp_, cutemp_, memsize, cudaMemcpyDefault, fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaMemcpyAsync failed with error code " << success << " (" << cudaGetErrorString((cudaError_t)success) << ")" << endl;
            parallel.abortForce();
          }

          success = cudaStreamSynchronize(fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
            parallel.abortForce();
          }
        }

        if (comp + 1 < components_ || alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : Transpose");
          b_transpose_back_0_1(cutemp1_, cutemp_,r2cSize_,r2cSizeLocal_as_,rSizeLocal_[2],rSizeLocal_[1],parallel.grid_size()[1], fft_stream);
          b_implement_0(&cutemp1_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]], cutemp_,r2cSize_,rSizeLocal_[2],rSizeLocal_[1], fft_stream);
          nvtxRangePop();

          nvtxRangePushA("Latfield2d::PlanFFT::execute : C2R transform");
          for(int l = 0;l< rSizeLocal_[2] ;l++)
          {
            p_cufft = &cutemp_[ l*r2cSize_*rSizeLocal_[1] ];
            p_out = &rData_[l*rJump_[2]*components_ + comp + (alignment_even ? 0 : 1)];

  #ifdef SINGLE
            cufft_status = cufftExecC2R(cubPlan_i_, p_cufft, (cufftReal*) p_out);
  #else
            cufft_status = cufftExecZ2D(cubPlan_i_, p_cufft, (cufftDoubleReal*) p_out);
  #endif
            if (cufft_status != CUFFT_SUCCESS)
            {
  #ifdef SINGLE
              cerr << "Latfield2d::PlanFFT::execute : cufftExecC2R failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #else
              cerr << "Latfield2d::PlanFFT::execute : cufftExecZ2D failed on component " << comp+(alignment_even ? 1 : 2) << "/" << components_ << " with error code " << cufft_status << " (" << cufftGetErrorString(cufft_status) << ")" << endl;
  #endif
              parallel.abortForce();
            }
          }
          nvtxRangePop();
        }

        if (comp + 1 < components_ || !alignment_even)
        {
          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved MPI AlltoAll and Scatter");
          MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
          MPI_Scatter(&temp_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp1_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
          MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
          nvtxRangePop();

          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved transpose");
          b_transpose_back_0_1(temp1_, temp_,r2cSize_,r2cSizeLocal_as_,rSizeLocal_[2],rSizeLocal_[1],parallel.grid_size()[1]);
          b_implement_0(&temp1_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]], temp_,r2cSize_,rSizeLocal_[2],rSizeLocal_[1]);
          nvtxRangePop();

          nvtxRangePushA("Latfield2d::PlanFFT::execute : interleaved C2R transform");
#pragma omp parallel for private(p_in, p_out)
          for(int l = 0;l< rSizeLocal_[2] ;l++)
          {
            p_in = &temp_[ l*r2cSize_*rSizeLocal_[1] ];
            p_out = &rData_[l*rJump_[2]*components_ + comp + (alignment_even ? 1 : 0)];
#ifdef SINGLE
            fftwf_execute_dft_c2r(bPlan_i_,p_in,p_out);
#else
            fftw_execute_dft_c2r(bPlan_i_,p_in,p_out);
#endif
          }
          nvtxRangePop();
        }

        if (comp + 2 >= components_ && (comp + 1 < components_ || alignment_even)) // defer synchronisation
        {
          success = cudaStreamSynchronize(fft_stream);
          if (success != cudaSuccess)
          {
            cerr << "Latfield2d::PlanFFT::execute : cudaStreamSynchronize failed" << endl;
            parallel.abortForce();
          }
        }
      }
    }

  }
  if(type_ == C2C)
  {

    		if(fft_type == FFT_FORWARD)
  			{

  			   //int i,j,k;
  				 //int comp;
  				 //int comm_rank;

  #ifdef SINGLE
  				 fftwf_complex *p_in;
  				 fftwf_complex *p_out;
  #else
  				 fftw_complex *p_in;
  				 fftw_complex *p_out;
  #endif

  				 for(int comp=0; comp<components_; comp++)
  				 {

  					for(int l = 0;l< rSizeLocal_[2] ;l++)
  					{

  					    p_in =  &cData_[rJump_[2]*l*components_ + comp];
  					    p_out = &temp_[l*rSizeLocal_[1]];
  #ifdef SINGLE
  					    fftwf_execute_dft(fPlan_i_,p_in,p_out);

  #else
  					    fftw_execute_dft(fPlan_i_,p_in,p_out);
  #endif
  					}


  				   MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);

  				   for(int i=0;i<parallel.grid_size()[1];i++)transpose_0_2(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[1]);

  #ifdef SINGLE
  				   fftwf_execute(fPlan_j_);
  #else
  				   fftw_execute(fPlan_j_);
  #endif


  				   MPI_Alltoall(temp_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1],MPI_DATA_PREC,temp1_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);

  				   for(int i=0;i<parallel.grid_size()[0];i++)transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]], rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[2]);


  				  for(int l = 0;l< rSizeLocal_[2] ;l++)
            //for(int l = 0;l< 13 ;l++)
  					{

  					    p_in  = &temp_[l*rSizeLocal_[1]];
  					    p_out = &kData_[rJump_[2]*l*components_ + comp];
                //p_out = &kData_[0];
  #ifdef SINGLE
  					    fftwf_execute_dft(fPlan_k_,p_in,p_out);
  #else
  					    fftw_execute_dft(fPlan_k_,p_in,p_out);
  #endif
  					}

            //p_in = NULL;
            //p_out = NULL;
  				 }

  			}
  			if(fft_type == FFT_BACKWARD)
  			{

  			   //int i,j,k;
  				 //int comp;
  				 //int comm_rank;

  #ifdef SINGLE
  				 fftwf_complex *p_in;
  				 fftwf_complex *p_out;
  #else
  				 fftw_complex *p_in;
  				 fftw_complex *p_out;
  #endif

  				 //STEP 1 : SAME AS STEP ONE OF FORWARD

  				for(int comp=0; comp<components_; comp++)
  				 {
  					for(int l = 0;l< rSizeLocal_[2] ;l++)
  					{

  					    p_in =  &kData_[rJump_[2]*l*components_ + comp];
  					    p_out = &temp_[l*rSizeLocal_[1]];
  #ifdef SINGLE
  					    fftwf_execute_dft(bPlan_k_,p_in,p_out);

  #else
  					    fftw_execute_dft(bPlan_k_,p_in,p_out);
  #endif
  					}


  				  //step 2 : same as step 4 of forward


  				  //MPI_Alltoall(temp_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1],MPI_DATA_PREC,temp1_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);

  				   //for(i=0;i<parallel.grid_size()[0];i++)transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]], rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[2]);


  #ifdef SINGLE
  				   fftwf_execute(bPlan_j_);
  #else
  				   fftw_execute(bPlan_j_);
  #endif


  				   MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);

  				   for(int i=0;i<parallel.grid_size()[1];i++)transpose_0_2(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[1]);


  				    for(int l = 0;l< rSizeLocal_[2] ;l++)
  					{

  					    p_in = &temp_[l*rSizeLocal_[1]];
  					    p_out =  &cData_[rJump_[2]*l*components_ + comp];

  #ifdef SINGLE
  					    fftwf_execute_dft(bPlan_i_,p_in,p_out);
                //fftwf_execute(bPlan_i_);
  #else
  					    fftw_execute_dft(bPlan_i_,p_in,p_out);
  #endif
  					}
  				 }

  			}
      }

  success = cudaStreamDestroy(fft_stream);
  if (success != cudaSuccess)
  {
    cerr << "Latfield2d::PlanFFT::execute : cudaStreamDestroy failed" << endl;
    parallel.abortForce(); 
  }
}

//transposition function
#ifdef SINGLE

/////
template<class compType>
void PlanFFT<compType>::transpose_0_2( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k)
{
  //int i,j,k;
  #pragma omp parallel for collapse(3)
  for(int i=0;i<dim_i;i++)
  {
    for(int j=0;j<dim_j;j++)
    {
      for(int k=0;k<dim_k;k++)
      {

        out[k+dim_k*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
        out[k+dim_k*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];

      }
    }
  }

}

__global__ void transpose_0_2_kernel(cufftComplex * in, cufftComplex * out, int dim_i, int dim_j, int dim_k)
{
  int i = blockIdx.x;
  int j = blockIdx.y;

  for (int k = threadIdx.x; k < dim_k; k += 128)
  {
    out[k + dim_k * (j + i * dim_j)] = in[i + dim_i * (j + k * dim_j)];
  }
}

template<class compType>
void PlanFFT<compType>::transpose_0_2( cufftComplex * in, cufftComplex * out, int dim_i, int dim_j, int dim_k, cudaStream_t &stream)
{
  transpose_0_2_kernel<<<dim3(dim_i, dim_j), 128, 0, stream>>>(in, out, dim_i, dim_j, dim_k);
}

template<class compType>
void PlanFFT<compType>::transpose_0_2_last_proc( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k)
{
  //int i,j,k;
  #pragma omp parallel for collapse(3)
  for(int i=0;i<dim_i;i++)
  {
    for(int j=0;j<dim_j;j++)
    {
      for(int k=0;k<dim_k;k++)
      {
        out[k+(dim_k+1)*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
        out[k+(dim_k+1)*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];
      }
    }
  }
}

__global__ void transpose_0_2_last_proc_kernel(cufftComplex * in, cufftComplex * out, int dim_i, int dim_j, int dim_k)
{
  int i = blockIdx.x;
  int j = blockIdx.y;

  for (int k = threadIdx.x; k < dim_k; k += 128)
  {
    out[k + (dim_k + 1) * (j + i * dim_j)] = in[i + dim_i * (j + k * dim_j)];
  }
}

template<class compType>
void PlanFFT<compType>::transpose_0_2_last_proc( cufftComplex * in, cufftComplex * out,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream)
{
  transpose_0_2_last_proc_kernel<<<dim3(dim_i, dim_j), 128, 0, stream>>>(in, out, dim_i, dim_j, dim_k);
}

template<class compType>
void PlanFFT<compType>::implement_local_0_last_proc( fftwf_complex * in, fftwf_complex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size)
{
  //int i_in,i_out,j,rank;
  #pragma omp parallel for collapse(3)
  for(int i_in=0;i_in<proc_dim_i;i_in++)
  {
    for(int j=0;j<proc_dim_j;j++)
    {
      for(int rank=0;rank<proc_size;rank++)
      {
        int i_out=i_in+rank*proc_dim_i;
        out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][0]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][0];
        out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][1]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][1];
      }
    }
  }
}

__global__ void implement_local_0_last_proc_kernel(cufftComplex * in, cufftComplex * out, int proc_dim_i, int proc_dim_j, int proc_dim_k)
{
  int rank = blockIdx.x;
  int j = blockIdx.y;

  for (int i_in = threadIdx.x; i_in < proc_dim_i; i_in += 128)
  {
    int i_out = i_in + rank * proc_dim_i;
    out[proc_dim_k + (proc_dim_k + 1) * (j + i_out * proc_dim_j)] = in[i_in + proc_dim_i * (j + proc_dim_j * rank)];
  }
}

template<class compType>
void PlanFFT<compType>::implement_local_0_last_proc( cufftComplex * in, cufftComplex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size, cudaStream_t &stream)
{
  implement_local_0_last_proc_kernel<<<dim3(proc_size, proc_dim_j), 128, 0, stream>>>(in, out, proc_dim_i, proc_dim_j, proc_dim_k);
}

template<class compType>
void PlanFFT<compType>::transpose_1_2(fftwf_complex * in , fftwf_complex * out ,int dim_i,int dim_j ,int dim_k )
{
  //int i,j,k;
  #pragma omp parallel for collapse(3)
  for(int i=0;i<dim_i;i++)
  {
    for(int j=0;j<dim_j;j++)
    {
      for(int k=0;k<dim_k;k++)
      {
        out[i+dim_i*(k+j*dim_k)][0]=in[i+dim_i*(j+k*dim_j)][0];
        out[i+dim_i*(k+j*dim_k)][1]=in[i+dim_i*(j+k*dim_j)][1];
      }
    }
  }
}

__global__ void transpose_1_2_kernel(cufftComplex * in, cufftComplex * out, int dim_i, int dim_j, int dim_k)
{
  int i = blockIdx.x;
  int j = blockIdx.y;

  for (int k = threadIdx.x; k < dim_k; k += 128)
  {
    out[i + dim_i * (k + j * dim_k)] = in[i + dim_i * (j + k * dim_j)];
  }
}

template<class compType>
void PlanFFT<compType>::transpose_1_2(cufftComplex * in , cufftComplex * out  ,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream)
{
  transpose_1_2_kernel<<<dim3(dim_i, dim_j), 128, 0, stream>>>(in, out, dim_i, dim_j, dim_k);
}

template<class compType>
void PlanFFT<compType>::transpose_back_0_3( fftwf_complex * in, fftwf_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components, int comp)
{
  //int i,j,k,l, i_t, j_t, k_t;
  int r2c_halo = r2c + 2*halo;
  int local_size_k_halo = local_size_k + 2*halo;

  #pragma omp parallel for collapse(4) default(shared)
  for (int i=0;i<local_r2c;i++)
  {
    for(int k=0;k<local_size_k;k++)
    {
      for(int j=0;j<local_size_j;j++)
      {
        for(int l=0;l<proc_size;l++)
        {
          int i_t = i + l*local_r2c;
          int j_t = j ;
          int k_t = k ;
          out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
          out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
        }
      }
    }
  }
}

template<class compType>
void PlanFFT<compType>::implement_0(fftwf_complex * in, fftwf_complex * out,int r2c_size,int local_size_j,int local_size_k, int halo,int components, int comp)
{
  //int i,j,k;
  int i=r2c_size-1;
  int r2c_halo = r2c_size + 2*halo;
  int local_size_k_halo = local_size_k + 2*halo;

  #pragma omp parallel for collapse(2) default(shared)
  for(int j=0;j<local_size_j;j++)
  {
    for(int k=0;k<local_size_k;k++)
    {
      out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][0]=in[j + local_size_j *k][0];
      out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][1]=in[j + local_size_j *k][1];
    }
  }

}


template<class compType>
void PlanFFT<compType>::b_arrange_data_0(fftwf_complex *in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k, int khalo, int components, int comp)
{
  //int i,j,k;
  int jump_i=(dim_i+ 2 *khalo);
  int jump_j=dim_j+ 2 *khalo;

  #pragma omp parallel for collapse(3) default(shared)
  for(int i=0;i<dim_i;i++)
  {
    for(int j=0;j<dim_j;j++)
    {
      for(int k=0;k<dim_k;k++)
      {
        out[j + dim_j * (k + dim_k * i)][0]=in[comp+components*(i + jump_i * (j + jump_j*k))][0];
        out[j + dim_j * (k + dim_k * i)][1]=in[comp+components*(i + jump_i * (j + jump_j*k))][1];
      }
    }
  }

}

template<class compType>
void PlanFFT<compType>::b_transpose_back_0_1( fftwf_complex * in, fftwf_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size)
{
  //int i,j,k,l, i_t, j_t, k_t;

  #pragma omp parallel for collapse(4) default(shared)
  for (int i=0;i<local_r2c;i++)
  {
    for(int k=0;k<local_size_k;k++)
    {
      for(int j=0;j<local_size_j;j++)
      {
        for(int l=0;l<proc_size;l++)
        {
          int i_t = i + l*local_r2c;
          int j_t = j ;
          int k_t = k ;
          out[i_t + r2c * (k_t + local_size_k * j_t)][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
          out[i_t + r2c * (k_t + local_size_k * j_t)][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
        }
      }
    }
  }
}

__global__ void b_transpose_back_0_1_kernel(cufftComplex * in, cufftComplex * out, int r2c, int local_r2c, int local_size_j, int local_size_k, int proc_size)
{
  int j = blockIdx.x;
  int k = blockIdx.y;

  //for (int l = 0; l < proc_size; l++)
  for (int i_t = threadIdx.x; i_t < local_r2c * proc_size; i_t += 128)
  {
    int i = i_t % local_r2c;
    int l = i_t / local_r2c;
    out[i_t + r2c * (k + local_size_k * j)] = in[i + local_r2c * (j + local_size_j * (k + local_size_k * l))];
  }
}

template<class compType>
void PlanFFT<compType>::b_transpose_back_0_1( cufftComplex * in, cufftComplex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size, cudaStream_t &stream)
{
  b_transpose_back_0_1_kernel<<<dim3(local_size_j, local_size_k), 128, 0, stream>>>(in, out, r2c, local_r2c, local_size_j, local_size_k, proc_size);
}

template<class compType>
void PlanFFT<compType>::b_implement_0(fftwf_complex * in, fftwf_complex * out,int r2c_size,int local_size_j,int local_size_k)
{
  //int i,j,k;
  int i=r2c_size-1;

  #pragma omp parallel for collapse(2) default(shared)
  for(int j=0;j<local_size_j;j++)
  {
    for(int k=0;k<local_size_k;k++)
    {
      out[i + r2c_size * (k + local_size_k *j)][0]=in[j + local_size_j *k][0];
      out[i + r2c_size * (k + local_size_k *j)][1]=in[j + local_size_j *k][1];
    }
  }
}

__global__ void b_implement_0_kernel(cufftComplex * in, cufftComplex * out, int r2c_size, int local_size_j, int local_size_k)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx < local_size_j * local_size_k)
  {
    int j = idx % local_size_j;
    int k = idx / local_size_j;
    out[r2c_size - 1 + r2c_size * (k + local_size_k * j)] = in[idx];
  }
}

template<class compType>
void PlanFFT<compType>::b_implement_0(cufftComplex * in, cufftComplex * out, int r2c_size, int local_size_j, int local_size_k, cudaStream_t &stream)
{
  b_implement_0_kernel<<<(local_size_j * local_size_k + 127) / 128, 128, 0, stream>>>(in, out, r2c_size, local_size_j, local_size_k);
}

#endif

#ifndef SINGLE

/////
template<class compType>
void PlanFFT<compType>::transpose_0_2( fftw_complex * in, fftw_complex * out,int dim_i,int dim_j ,int dim_k)
{
  //int i,j,k;

  #pragma omp parallel for collapse(3)
  for(int i=0;i<dim_i;i++)
  {
    for(int j=0;j<dim_j;j++)
    {
      for(int k=0;k<dim_k;k++)
      {

        out[k+dim_k*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
        out[k+dim_k*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];

      }
    }
  }

}

__global__ void transpose_0_2_kernel(cufftDoubleComplex * in, cufftDoubleComplex * out, int dim_i, int dim_j, int dim_k)
{
  int i = blockIdx.x;
  int j = blockIdx.y;

  for (int k = threadIdx.x; k < dim_k; k += 128)
  {
    out[k + dim_k * (j + i * dim_j)] = in[i + dim_i * (j + k * dim_j)];
  }
}

template<class compType>
void PlanFFT<compType>::transpose_0_2( cufftDoubleComplex * in, cufftDoubleComplex * out, int dim_i, int dim_j, int dim_k, cudaStream_t &stream)
{
  transpose_0_2_kernel<<<dim3(dim_i, dim_j), 128, 0, stream>>>(in, out, dim_i, dim_j, dim_k);
}

template<class compType>
void PlanFFT<compType>::transpose_0_2_last_proc( fftw_complex * in, fftw_complex * out,int dim_i,int dim_j ,int dim_k)
{
  //int i,j,k;

  #pragma omp parallel for collapse(3)
  for(int i=0;i<dim_i;i++)
  {
    for(int j=0;j<dim_j;j++)
    {
      for(int k=0;k<dim_k;k++)
      {
        out[k+(dim_k+1)*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
        out[k+(dim_k+1)*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];
      }
    }
  }
}

__global__ void transpose_0_2_last_proc_kernel(cufftDoubleComplex * in, cufftDoubleComplex * out, int dim_i, int dim_j, int dim_k)
{
  int i = blockIdx.x;
  int j = blockIdx.y;

  for (int k = threadIdx.x; k < dim_k; k += 128)
  {
    out[k + (dim_k + 1) * (j + i * dim_j)] = in[i + dim_i * (j + k * dim_j)];
  }
}

template<class compType>
void PlanFFT<compType>::transpose_0_2_last_proc( cufftDoubleComplex * in, cufftDoubleComplex * out,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream)
{
  transpose_0_2_last_proc_kernel<<<dim3(dim_i, dim_j), 128, 0, stream>>>(in, out, dim_i, dim_j, dim_k);
}

template<class compType>
void PlanFFT<compType>::implement_local_0_last_proc( fftw_complex * in, fftw_complex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size)
{
  //int i_in,i_out,j,rank;

  #pragma omp parallel for collapse(3)
  for(int i_in=0;i_in<proc_dim_i;i_in++)
  {
    for(int j=0;j<proc_dim_j;j++)
    {
      for(int rank=0;rank<proc_size;rank++)
      {
        int i_out=i_in+rank*proc_dim_i;
        out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][0]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][0];
        out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][1]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][1];
      }
    }
  }
}

__global__ void implement_local_0_last_proc_kernel(cufftDoubleComplex * in, cufftDoubleComplex * out, int proc_dim_i, int proc_dim_j, int proc_dim_k)
{
  int rank = blockIdx.x;
  int j = blockIdx.y;

  for (int i_in = threadIdx.x; i_in < proc_dim_i; i_in += 128)
  {
    int i_out = i_in + rank * proc_dim_i;
    out[proc_dim_k + (proc_dim_k + 1) * (j + i_out * proc_dim_j)] = in[i_in + proc_dim_i * (j + proc_dim_j * rank)];
  }
}

template<class compType>
void PlanFFT<compType>::implement_local_0_last_proc( cufftDoubleComplex * in, cufftDoubleComplex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size, cudaStream_t &stream)
{
  implement_local_0_last_proc_kernel<<<dim3(proc_size, proc_dim_j), 128, 0, stream>>>(in, out, proc_dim_i, proc_dim_j, proc_dim_k);
}

template<class compType>
void PlanFFT<compType>::transpose_1_2(fftw_complex * in , fftw_complex * out ,int dim_i,int dim_j ,int dim_k )
{
  //int i,j,k;

  #pragma omp parallel for collapse(3)
  for(int i=0;i<dim_i;i++)
  {
    for(int j=0;j<dim_j;j++)
    {
      for(int k=0;k<dim_k;k++)
      {
        out[i+dim_i*(k+j*dim_k)][0]=in[i+dim_i*(j+k*dim_j)][0];
        out[i+dim_i*(k+j*dim_k)][1]=in[i+dim_i*(j+k*dim_j)][1];
      }
    }
  }
}

__global__ void transpose_1_2_kernel(cufftDoubleComplex * in, cufftDoubleComplex * out, int dim_i, int dim_j, int dim_k)
{
  int i = blockIdx.x;
  int j = blockIdx.y;

  for (int k = threadIdx.x; k < dim_k; k += 128)
  {
    out[i + dim_i * (k + j * dim_k)] = in[i + dim_i * (j + k * dim_j)];
  }
}

template<class compType>
void PlanFFT<compType>::transpose_1_2(cufftDoubleComplex * in , cufftDoubleComplex * out  ,int dim_i,int dim_j ,int dim_k, cudaStream_t &stream)
{
  transpose_1_2_kernel<<<dim3(dim_i, dim_j), 128, 0, stream>>>(in, out, dim_i, dim_j, dim_k);
}

template<class compType>
void PlanFFT<compType>::transpose_back_0_3( fftw_complex * in, fftw_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components, int comp)
{
  //int i,j,k,l, i_t, j_t, k_t;
  int r2c_halo = r2c + 2*halo;
  int local_size_k_halo = local_size_k + 2*halo;

  #pragma omp parallel for collapse(4) default(shared)
  for (int i=0;i<local_r2c;i++)
  {
    for(int k=0;k<local_size_k;k++)
    {
      for(int j=0;j<local_size_j;j++)
      {
        for(int l=0;l<proc_size;l++)
        {
          int i_t = i + l*local_r2c;
          int j_t = j ;
          int k_t = k ;
          out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
          out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
        }
      }
    }
  }
}

template<class compType>
void PlanFFT<compType>::implement_0(fftw_complex * in, fftw_complex * out,int r2c_size,int local_size_j,int local_size_k, int halo,int components, int comp)
{
  // int i,j,k;
  int i=r2c_size-1;
  int r2c_halo = r2c_size + 2*halo;
  int local_size_k_halo = local_size_k + 2*halo;

  #pragma omp parallel for collapse(2) default(shared)
  for(int j=0;j<local_size_j;j++)
  {
    for(int k=0;k<local_size_k;k++)
    {
      out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][0]=in[j + local_size_j *k][0];
      out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][1]=in[j + local_size_j *k][1];
    }
  }

}


template<class compType>
void PlanFFT<compType>::b_arrange_data_0(fftw_complex *in, fftw_complex * out,int dim_i,int dim_j ,int dim_k, int khalo, int components, int comp)
{
  //int i,j,k;
  int jump_i=(dim_i+ 2 *khalo);
  int jump_j=dim_j+ 2 *khalo;

  #pragma omp parallel for collapse(3) default(shared)
  for(int i=0;i<dim_i;i++)
  {
    for(int j=0;j<dim_j;j++)
    {
      for(int k=0;k<dim_k;k++)
      {
        out[j + dim_j * (k + dim_k * i)][0]=in[comp+components*(i + jump_i * (j + jump_j*k))][0];
        out[j + dim_j * (k + dim_k * i)][1]=in[comp+components*(i + jump_i * (j + jump_j*k))][1];
      }
    }
  }

}

template<class compType>
void PlanFFT<compType>::b_transpose_back_0_1( fftw_complex * in, fftw_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size)
{
  //int i,j,k,l, i_t, j_t, k_t;

  #pragma omp parallel for collapse(4)
  for (int i=0;i<local_r2c;i++)
  {
    for(int k=0;k<local_size_k;k++)
    {
      for(int j=0;j<local_size_j;j++)
      {
        for(int l=0;l<proc_size;l++)
        {
          int i_t = i + l*local_r2c;
          int j_t = j ;
          int k_t = k ;
          out[i_t + r2c * (k_t + local_size_k * j_t)][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
          out[i_t + r2c * (k_t + local_size_k * j_t)][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
        }
      }
    }
  }
}

__global__ void b_transpose_back_0_1_kernel(cufftDoubleComplex * in, cufftDoubleComplex * out, int r2c, int local_r2c, int local_size_j, int local_size_k, int proc_size)
{
  int j = blockIdx.x;
  int k = blockIdx.y;

  //for (int l = 0; l < proc_size; l++)
  for (int i_t = threadIdx.x; i_t < local_r2c * proc_size; i_t += 128)
  {
    int i = i_t % local_r2c;
    int l = i_t / local_r2c;
    out[i_t + r2c * (k + local_size_k * j)] = in[i + local_r2c * (j + local_size_j * (k + local_size_k * l))];
  }
}

template<class compType>
void PlanFFT<compType>::b_transpose_back_0_1( cufftDoubleComplex * in, cufftDoubleComplex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size, cudaStream_t &stream)
{
  b_transpose_back_0_1_kernel<<<dim3(local_size_j, local_size_k), 128, 0, stream>>>(in, out, r2c, local_r2c, local_size_j, local_size_k, proc_size);
}

template<class compType>
void PlanFFT<compType>::b_implement_0(fftw_complex * in, fftw_complex * out,int r2c_size,int local_size_j,int local_size_k)
{
  //int i,j,k;
  int i=r2c_size-1;

  #pragma omp parallel for collapse(2) default(shared)
  for(int j=0;j<local_size_j;j++)
  {
    for(int k=0;k<local_size_k;k++)
    {
      out[i + r2c_size * (k + local_size_k *j)][0]=in[j + local_size_j *k][0];
      out[i + r2c_size * (k + local_size_k *j)][1]=in[j + local_size_j *k][1];
    }
  }

}

__global__ void b_implement_0_kernel(cufftDoubleComplex * in, cufftDoubleComplex * out, int r2c_size, int local_size_j, int local_size_k)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx < local_size_j * local_size_k)
  {
    int j = idx % local_size_j;
    int k = idx / local_size_j;
    out[r2c_size - 1 + r2c_size * (k + local_size_k * j)] = in[idx];
  }
}

template<class compType>
void PlanFFT<compType>::b_implement_0(cufftDoubleComplex * in, cufftDoubleComplex * out, int r2c_size, int local_size_j, int local_size_k, cudaStream_t &stream)
{
  b_implement_0_kernel<<<(local_size_j * local_size_k + 127) / 128, 128, 0, stream>>>(in, out, r2c_size, local_size_j, local_size_k);
}

#endif


#endif


#endif
