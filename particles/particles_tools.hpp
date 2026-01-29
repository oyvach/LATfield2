#ifndef PARTICLES_TOOLS_HPP
#define PARTICLES_TOOLS_HPP


#ifndef DOXYGEN_SHOULD_SKIP_THIS
LATfield2::Real get_lattice_resolution(int npts[3],LATfield2::Real boxSize[3])
{
  LATfield2::Real latRes[3];

  for(int i =0;i<3;i++)latRes[i]=boxSize[i]/ (LATfield2::Real)npts[i];

  if(latRes[0]==latRes[1] && latRes[0]==latRes[2])return latRes[0];
  else{
    COUT<< "wrong physical box size and lattice size, relosution must be same in each dimensions"<<endl;
    COUT << " latRes[0]="<<latRes[0]<<" latRes[1]="<<latRes[1]<<" latRes[2]="<<latRes[2]<<endl;
    COUT << " boxSize[0]="<<boxSize[0]<<" boxSize[1]="<<boxSize[1]<<" boxSize[2]="<<boxSize[2]<<endl;
    COUT << " npts[0]="<<npts[0]<<" npts[1]="<<npts[1]<<" npts[2]="<<npts[2]<<endl;
    COUT<< " Exiting... "<<endl;
    exit(222);
      return -1;
  }
}


#define CREATE_MEMBER_DETECTOR(X)            \
  template<typename T> struct has_##X {      \
    struct Fallback {int X; };               \
    struct Derived : T, Fallback { };        \
    template<typename C, C> struct ChT;      \
    template<typename C> static char (&f(ChT<int Fallback::*, &C::X>*))[1]; \
    template<typename C> static char (&f(...))[2]; \
    static bool const value = sizeof(f<Derived>(0)) == 2; \
  };


#define CREATE_MEMBER_DETECTOR_MAXI(X)          \
    template<typename T> struct has_maxi_##X {      \
    struct Fallback {int X; };               \
    struct Derived : T, Fallback { };        \
    template<typename C, C> struct ChT;      \
    template<typename C,typename CC> int (f(ChT<int Fallback::*, &C::X>*)){return -1;}; \
    template<typename C,typename CC> int (f(...)){return offsetof(T,X);}; \
    int gos(){return f<Derived,T>(0);} \
};

#endif

#endif