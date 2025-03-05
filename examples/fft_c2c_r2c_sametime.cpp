/*! file poissonSolver.cpp
    Created by David Daverio.

    A simple example of LATfield2d usage. This exec solve the poisson equation in fourier space.


 */
#include <vector>
#include <iostream>
#include "LATfield2.hpp"
using namespace LATfield2;

int main(int argc, char **argv)
{
  int n, m;
  int BoxSize = 16;
  int halo = 0;
  int khalo = 1;
  int dim = 3;
  int comp = 1;

  for (int i = 1; i < argc; i++)
  {
    if (argv[i][0] != '-')
      continue;
    switch (argv[i][1])
    {
    case 'n':
      n = atoi(argv[++i]);
      break;
    case 'm':
      m = atoi(argv[++i]);
      break;
    case 'b':
      BoxSize = atoi(argv[++i]);
      break;
    }
  }

  parallel.initialize(n, m);

  Lattice lat;
  lat.initialize(3, BoxSize, halo);

  // Real to complex fourier transform

  Lattice latKreal;
  latKreal.initializeRealFFT(lat, khalo);

  Site x(lat);
  rKSite kReal(latKreal);

  Field<Real> phi;
  phi.initialize(lat, comp);

  Field<Imag> phiK;
  phiK.initialize(latKreal, comp);

  PlanFFT<Imag> planReal(&phi, &phiK);

  Lattice latKcomplex;
  latKcomplex.initializeComplexFFT(lat, khalo);
  cKSite kComplex(latKcomplex);

  Field<Imag> rho;
  rho.initialize(lat, comp);

  Field<Imag> rhoK;
  rhoK.initialize(latKcomplex, comp);

  PlanFFT<Imag> planComplex(&rho, &rhoK);

  double line = (double)phi.lattice().size(0);
  double l3 = line * line * line;
  COUT << "Line = " << line << endl;

  for (x.first(); x.test(); x.next())
  {
    phi(x) = 1.0;
  }
  double avreal = 0., avcom = 0.;
  for (x.first(); x.test(); x.next())
  {
    avreal += phi(x);
  }
  parallel.sum<double>(avreal);
  avreal *= 1./l3;
  COUT << "Average real at start = " << avreal << endl;
  COUT << "Starting r2c fouriers" << endl;
  planReal.execute(FFT_FORWARD);
  planReal.execute(FFT_BACKWARD);
  avreal = 0.;
  for (x.first(); x.test(); x.next())
  {
    avreal += phi(x);
  }
  parallel.sum<double>(avreal);
  avreal *= 1./l3/l3;
  COUT << "Average real = " << avreal << endl;

  

  for (x.first(); x.test(); x.next())
  {
    rho(x).real() = 1.0;
    rho(x).imag() = 0.0;
  }

  planComplex.execute(FFT_FORWARD);
  planComplex.execute(FFT_BACKWARD);

  avreal = 0.;
  avcom = 0.;
  for (x.first(); x.test(); x.next())
  {
    avreal += phi(x);
    avcom += rho(x).real();
  }

  parallel.sum<double>(avreal);
  avreal *= 1./l3/l3;
  parallel.sum<double>(avcom);
  avcom *= 1./l3/line/line; // conclusion that c2c Fourier only multiplies line2, while r2c multiplies line3?
  COUT << "Average real = " << avreal << " Average complex = " << avcom << endl;
  planReal.execute(FFT_FORWARD);
  planReal.execute(FFT_BACKWARD);
  COUT << "Finished" << endl;


}
