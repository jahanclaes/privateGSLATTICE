//#define checkderivs
#include <iostream>
#include <fstream>
#include <cmath>
#include "SystemClass.h"
#include "WaveFunction.h" 
#include "Jastrow.h"
void JastrowClass::MakeUniformState()
{
  for (int i=0;i<vij.size();i++){
    vij(i)=0.0;
  }
}

void JastrowClass::Init(SystemClass &system)
{
  //Read the jastrow parameters
  GetParameters();

  Ti.resize(system.x.size());
  InitTi(system);

#ifdef debug

  int i,j,ix,iy,jx,jy;

  for(i = 0; i < system.x.size(); i++)
  {
    ix = i % system.N1; iy = i / system.N1;

    for(j = 0; j < system.x.size(); j++)
    {
      jx = j % system.N1; jy = j / system.N1;
      cout << "(" << ix << " " << iy << ") " << "(" << jx << " " << jy << ") "; 
      cout << GetReducedIndex(system,i,j) << endl;
    }
  }

assert(1 == 2);

#endif
  cout << "Complete th initialization of Jastrow type" << endl;
}

complex<double> JastrowClass::logevaluate(SystemClass &system,int &sign)
{}

void JastrowClass::CheckDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start, int stop)
{
  int i,j,k,wj,wk;
  double var_tmp,wij,PJ_new,PJ_old;
  
  for(i = start; i < stop; i++)
  {
    wij = 0.;

    for(j = 0; j < system.x.size(); j++)
    {
      wj = abs(system.x(j)) - 1;

      for(k = 0; k < system.x.size(); k++)
      {
        wk   = abs(system.x(k)) - 1;
        wij += vij(GetReducedIndex(system,j,k)) * (double)(wj * wk);
      }
    }

    PJ_old = exp(0.5 * wij);

    var_tmp = vij(reduced_index(i-start));
    vij(reduced_index(i-start)) += 0.0001;

    wij = 0.;

    for(j = 0; j < system.x.size(); j++)
    {
      wj = abs(system.x(j)) - 1;

      for(k = 0; k < system.x.size(); k++)
      {
        wk   = abs(system.x(k)) - 1;
        wij += vij(GetReducedIndex(system,j,k)) * (double)(wj * wk);
      }
    }

    PJ_new = exp(0.5 * wij);

    derivs(i) = (PJ_new - PJ_old) * 10000.0 / PJ_old;

    vij(reduced_index(i-start)) = var_tmp;
  }
}

void JastrowClass::AllDerivs(SystemClass &system, Array<complex<double>,1> 
                                          &derivs, int start, int stop)
{
  int i,j,k,rij,wj,wk;

  for(i = start; i < stop; i++)
  {
    rij = reduced_index(i-start);
    derivs(i).real(0);
    //    derivs(i).real() = 0.;

    for(j = 0; j < system.x.size(); j++)
    {
      wj = abs(system.x(j)) - 1;

      for(k = 0; k < system.x.size(); k++)
      {
        if(GetReducedIndex(system,j,k) == rij)
        {
          wk = abs(system.x(k)) - 1;
	  //    derivs(i).real() += (double) (wj * wk);
	  derivs(i).real((double) (wj * wk)+derivs(i).real());
        }
      }
    }

    derivs(i) *= 0.5;
  }

#ifdef checkderivs

  Array<complex<double>,1> derivs_tmp;
  derivs_tmp.resize(derivs.size());
  CheckDerivs(system,derivs_tmp,start,stop);
  
  for(i = start; i < stop; i++)
    cout << "Check derivs " << derivs(i) << " " << derivs_tmp(i) << endl;

  derivs_tmp.free();

#endif
 
}

void JastrowClass::SetParam(int i, double param)
{
  vij(reduced_index(i)) = param;
}

double JastrowClass::GetParam(int i)
{
  return vij(reduced_index(i));
}

double JastrowClass::GetParam_real(int i)
{
  return GetParam(i);
}

void JastrowClass::SetParam_real(int i, double param)
{
  SetParam(i,param);
}
double JastrowClass::GetParam_imag(int i)
{
  return 0.;
}
void JastrowClass::SetParam_imag(int i, double param){}

complex<double> JastrowClass::evaluate(SystemClass &system){}
complex<double> JastrowClass::evaluateRatio(SystemClass &system,int swap1, int swap2){}
complex<double> JastrowClass::evaluateRatio_check(SystemClass &system, int swap1, int swap2){}
void JastrowClass::UpdateDets(SystemClass &system,int swap1, int swap2){};

void JastrowClass::UpdateDets(SystemClass &system,int site, int end_site,int spin)
{
  //update determinant and update Ti
  //cout << "Update T matrix" << endl;
  //int i,s = abs(spin);
  int i;

  for(i = 0; i < system.x.size(); i++)
    //Ti(i) += (double)s * (vij(GetReducedIndex(system,i,end_site)) - vij(GetReducedIndex(system,i,site)));
    Ti(i) += (vij(GetReducedIndex(system,i,end_site)) - vij(GetReducedIndex(system,i,site)));
}

complex<double> JastrowClass::evaluateRatio_check(SystemClass &system, int start,int stop,int spin)
{
  int i,j,ni,nj;
  double wij = 0.,PJ_new,PJ_old;
  complex<double> ratio;

  for(j = 0; j < system.x.size(); j++)
  {
    nj = abs(system.x(j));

    for(i = 0; i < system.x.size(); i++)
    {
      ni = abs(system.x(i));
      wij += vij(GetReducedIndex(system,j,i)) * ni * nj;
    }
  }

  PJ_new = exp(0.5 * wij);

  system.Move(stop,start,spin);
  wij = 0.;

  for(j = 0; j < system.x.size(); j++)
  {
    nj = abs(system.x(j));
    for(i = 0; i < system.x.size(); i++)
    {
      ni = abs(system.x(i));
      wij += vij(GetReducedIndex(system,j,i)) * ni * nj;
    }
  }

  PJ_old = exp(0.5 * wij);
  
  ratio = PJ_new/PJ_old;

  system.Move(start,stop,spin);  

  return ratio;
 
}        

complex<double> JastrowClass::evaluateRatio(SystemClass &system,int start,int stop,int spin)
{
  //int s;
  double t_l,t_k,v_ll,v_lk;
  complex<double> ratio;

  //s   = abs(spin);
  t_l = Ti(stop); t_k = Ti(start);
  v_ll = vij(0);
  v_lk = vij(GetReducedIndex(system,start,stop));
  
  //ratio = exp((double)s * (t_l - t_k) + v_ll - v_lk);
  ratio = exp((t_l - t_k) + v_ll - v_lk);

#ifdef checkratio

  complex<double> ratio_check;

  ratio_check = evaluateRatio_check(system,start,stop,spin);

  cout << "Check ratio: " << ratio << " " << ratio_check << endl;

#endif

  return ratio;
 
}

void JastrowClass::Swap(int i, int j){}
void JastrowClass::Move(int site, int end_site, int s) {}
void JastrowClass::Reject(SystemClass &system,int site,int end_site,int spin) {}

// input a pair of indices and return an irreducible index
int JastrowClass::GetReducedIndex(SystemClass &system,int i,int j)
{
  int k,ix,iy,jx,jy,dx,dy,vx,vy,N1,N2,idx;

  N1 = system.N1; N2 = system.N2;
  ix = i % N1; iy = i / N1;
  jx = j % N1; jy = j / N1;

  dx = abs(ix-jx); dy = abs(iy-jy);

  vx  = ( (2 * dx) > N1 ) ? (N1 - dx) : dx;
  vy  = ( (2 * dy) > N2 ) ? (N2 - dy) : dy;
  idx = (vx > vy) ? vx + vy * N1 : vy + vx * N1;

  return idx;
}

void JastrowClass::InitTi(SystemClass &system)
{
  cout << "Init T matrix" << endl;
  int i,j;

  for(j = 0; j < system.x.size(); j++)
  {
    Ti(j) = 0.;

    for(i = 0; i < system.x.size(); i++)
    {
      //cout << j << " " << i << " " << GetReducedIndex(system,j,i) << endl;
      Ti(j) += vij(GetReducedIndex(system,j,i)) * 
               (double)(abs(system.x(i)) - 1);
    }
  }
}

//This is the initialize for the parameters
void JastrowClass::GetParameters()
{
  int number_reduced_index;
  ifstream infile;
  infile.open("jastrow.txt");
 
  infile >> NumParams;
  cout << "The number of jastrow parameters: " << NumParams << endl;

  reduced_index.resize(NumParams);
 
  for(int i = 0; i < NumParams; i++)
    infile >> reduced_index(i); 

  infile >> number_reduced_index;

 
  vij.resize(number_reduced_index);

  for(int i = 0; i < number_reduced_index; i++)
    infile >> vij(i); 

  infile.close();

  cout << "The Jastrow parameters:" << endl;
  for(int i = 0; i < NumParams; i++)
    cout << vij(reduced_index(i)) << " "; 
  cout << endl;

}

