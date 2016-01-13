//#define checkderivs
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

#include "SystemClass.h"
#include "WaveFunction.h" 
#include "BackFlow.h"
#include "SmartEigen.h"

void BackFlowClass::Init(SystemClass &system)
{
  // read the orbital from file
  eigs.resize(system.kList.size(),system.rList.size());
  GetEigs(system);

  //Initialize the spin config.
  spin.resize(2 * system.x.size());
  InitSpin(system);

  //Initialize the non-backflow retegular M
  Ma.resize(n_a,system.x.size());
  Mb.resize(n_b,system.x.size());
  InitM(system);

  eigs.free();

  //Initialize the square Da, Db, DA, DB 
  Da.Init(n_a); Db.Init(n_b);
  DA.Init(n_a); DB.Init(n_b);

  //Read the parameters \eta and \epsilon and neighbors
  GetParameters(system);

  //Initialize the backflow orbitals with the values of non-backflow.
  if(level > 0)
  {
    Ma_bf.resize(n_a,system.x.size());
    Mb_bf.resize(n_b,system.x.size());
  }
}

complex<double> BackFlowClass::logevaluate(SystemClass &system,int &sign)
{}

void BackFlowClass::CheckDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start, int stop)
{
  int i,j,k;
  complex<double> Det_a_new1,Det_b_new1,var_tmp;
  complex<double> Det_a_new2,Det_b_new2;
  complex<double> Det_a_old,Det_b_old;

  ConvertSpin(system);
  BackFlowCoord(system);
  FillDet(system,Da,Db,Ma_bf,Mb_bf);
  Det_a_old = Da.Det(); Det_b_old = Db.Det();
  
  for(i = start; i < stop; i++)
  {
    var_tmp = var(i);
    var(i).real() += 0.0001;

    BackFlowCoord(system);
    FillDet(system,Da,Db,Ma_bf,Mb_bf);
    Det_a_new1 = Da.Det(); Det_b_new1 = Db.Det();

    var(i).real() -= 0.0002;

    BackFlowCoord(system);
    FillDet(system,Da,Db,Ma_bf,Mb_bf);
    Det_a_new2 = Da.Det(); Det_b_new2 = Db.Det();

    derivs(i) = ( Det_a_new1 * Det_b_new1 - Det_a_new2 * Det_b_new2 ) * 5000.0;
    derivs(i) = derivs(i) / (Det_a_old * Det_b_old);

    var(i) = var_tmp;
  }
}

void BackFlowClass::AllDerivs(SystemClass &system, Array<complex<double>,1> 
                                          &derivs, int start, int stop)
{
  int i,j,k;
  complex<double> Det_a,Det_b,tmp,var_tmp;

  Array <complex<double>,2> Da_tmp,Db_tmp;

  ConvertSpin(system);

  DA.CalcAndSaveInverse();
  DB.CalcAndSaveInverse();

  Da_tmp.resize(n_a,n_a); Db_tmp.resize(n_b,n_b);

  //Compute d\det(\phi)/d\var/det(\phi)
  for(i = start; i < stop; i++)
  {
    derivs(i) = (0.,0.);
    var_tmp = var(i);
    
    var(i).real() = 1.; var(i).imag() = 0.;
    BackFlowCoord(system);
    FillDet(system,Da_tmp,Db_tmp,Ma_bf,Mb_bf);


    var(i).real() = 0.; var(i).imag() = 0.; 
    BackFlowCoord(system);
    FillDet(system,Da,Db,Ma_bf,Mb_bf);

    tmp = (0.,0.); 
    for(k = 0; k < Da_tmp.extent(0); k++)
      for(j = 0; j < Da_tmp.extent(1); j++)
        tmp += DA.MInverse(j,k) * (Da_tmp(k,j) - Da.M(k,j));

    derivs(i) = tmp;

    tmp = (0.,0.);    
    for(k = 0; k < Db_tmp.extent(0); k++)
      for(j = 0; j < Db_tmp.extent(1); j++)
        tmp += DB.MInverse(j,k) * (Db_tmp(k,j) - Db.M(k,j));

    derivs(i) += tmp;

    var(i) = var_tmp;
  }

#ifdef checkderivs

  if(level > 0)
  {
    Array<complex<double>,1> derivs_tmp;
    derivs_tmp.resize(derivs.size());
    CheckDerivs(system,derivs_tmp,start,stop);
  
    for(i = start; i < stop; i++)
      cout << "Check derivs " << derivs(i) << " " << derivs_tmp(i) << endl;

    derivs_tmp.free();
  }
  else
    cout << "No check for backflow derivative" << endl;

#endif
 
  Da_tmp.free(); Db_tmp.free(); 
}

void BackFlowClass::SetParam_real(int i, double param)
{
  var(i).real()=param;
}

void BackFlowClass::SetParam_imag(int i, double param)
{
  var(i).imag()=param;
}

double BackFlowClass::GetParam_real(int i)
{
  return var(i).real();
}

double BackFlowClass::GetParam_imag(int i)
{
  return var(i).imag();
}

complex<double> BackFlowClass::evaluateRatio(SystemClass &system,int start, int stop, int spin)
{
  complex <double> Det_a,Det_b,Det_a_new,Det_b_new,ratio;

  if(level > 0)
  {
    ConvertSpin(system);
    BackFlowCoord(system);
    FillDet(system,Da,Db,Ma_bf,Mb_bf);
  }
  else
    FillDet(system,Da,Db,Ma,Mb);

  DET_a_new = Da.Det(); DET_b_new = Db.Det();

#ifdef debug
  cout << "The new spin-up determinant: " << Det_a_new << endl;    
  cout << "The new spin-down determinant: " << Det_b_new << endl;
#endif

  ratio = (DET_a_new * DET_b_new) / (DET_a * DET_b);

#ifdef debug
  cout << "The ratio: " << ratio << endl;
#endif

  return ratio;
}

void BackFlowClass::CopyMatrix(Array<complex<double>,2> &Src,Array<complex<double>,2> &Dest)
{
  int i,k;

  #pragma omp parallel for private(i,k)
  for(k = 0; k < Dest.extent(0); k++)
    for(i = 0; i < Dest.extent(1); i++)
      Dest(k,i) = Src(k,i); 
}

void BackFlowClass::CopyMatrix(SmartEigen &Src,Array<complex<double>,2> &Dest)
{
  int i,k;

  #pragma omp parallel for private(i,k)
  for(k = 0; k < Dest.extent(0); k++)
    for(i = 0; i < Dest.extent(1); i++)
      Dest(k,i) = Src.M(k,i);
}

void BackFlowClass::CopyMatrix(SmartEigen &Src,SmartEigen &Dest,int size)
{
  int i,k;

  #pragma omp parallel for private(i,k)
  for(k = 0; k < size; k++)
    for(i = 0; i < size; i++)
      Dest.M(k,i) = Src.M(k,i);
}

void BackFlowClass::UpdateDets(SystemClass &system,int site, int end_site,int spin)
{
  int i,j;

  CopyMatrix(Da,DA,n_a);
  CopyMatrix(Db,DB,n_b);

  DET_a = DET_a_new; DET_b = DET_b_new;
}

void BackFlowClass::Reject(SystemClass &system,int site,int end_site,int spin)
{
}

void BackFlowClass::GetEigs(SystemClass &system)
{
  int i,j;
  ifstream infile;

  infile.open("eigs.txt");

  for (i = 0;i < eigs.extent(0); i++)
  {
    for (j = 0;j < eigs.extent(1);j++)
    {      
      infile >> eigs(i,j).real();
      //eigs(i,j).imag() = 0.;
      infile >> eigs(i,j).imag();
    }
  }

  infile.close();
  cout <<"Get orbitals successfully\n";

}

void BackFlowClass::InitSpin(SystemClass &system)
{
  int i;

  for(i = 0; i < system.x.size(); i++)
  {
    switch(system.x(i))
    {
      case 1:
      n_a++;break;

      case 2:
      n_a++; n_b++; break;

      case 0:
      break;

      case -1:
      n_b++; break;
    }
  }
  n_p = n_a + n_b;
  cout << "The number of up and down electrons:" 
       << n_a << " "<< n_b << " " << n_p << "\n";
}

void BackFlowClass::ConvertSpin(SystemClass &system)
{
  int i;

  for(i = 0; i < system.x.size(); i++)
  {
    switch(system.x(i))
    {
      case 1:
      spin(2 * i) = 1; spin(2 * i + 1) = 0;
      break;

      case 2:
      spin(2 * i) = 1;spin(2 * i + 1) = 1;
      break;

      case 0:
      spin(2 * i) = 0; spin(2 * i + 1) = 0;
      break;

      case -1:
      spin(2 * i) = 0; spin(2 * i + 1) = 1;
      break;
    }
  }
}

void BackFlowClass::FillDet(SystemClass &system,SmartEigen &D_a,SmartEigen &D_b,
Array<complex<double>,2> &M_a,Array<complex<double>,2> &M_b)
{
  int i,k,ja = 0,jb = 0;
 
  for(i = 0; i < system.x.size(); i++)
  {
    switch(system.x(i))
    {
      case 1:
        for(k = 0; k < M_a.extent(0); k++)  // loop in sites
          D_a.M(k,ja) = M_a(k,i);
        ja++; break;

      case 2:
        for(k = 0; k < M_a.extent(0); k++)  // loop in sites
          D_a.M(k,ja) = M_a(k,i);
     
        for(k = 0; k < M_b.extent(0); k++)
          D_b.M(k,jb) = M_b(k,i);
  
        ja++; jb++; break;

      case -1:
        for(k = 0; k < M_b.extent(0); k++)
          D_b.M(k,jb) = M_b(k,i);
   
        jb++; break;

      case 0:
        break;
    }
  }
}

void BackFlowClass::FillDet(SystemClass &system,Array<complex<double>,2> &D_a,
Array<complex<double>,2> &D_b, Array<complex<double>,2> &M_a,Array<complex<double>,2> &M_b)
{
  int i,k,ja = 0,jb = 0;

  for(i = 0; i < system.x.size(); i++)
  {
    switch(system.x(i))
    {
      case 1:
        for(k = 0; k < M_a.extent(0); k++)  // loop in sites
          D_a(k,ja) = M_a(k,i);
        ja++; break;

      case 2:
        for(k = 0; k < M_a.extent(0); k++)  // loop in sites
          D_a(k,ja) = M_a(k,i);

        for(k = 0; k < M_b.extent(0); k++)
          D_b(k,jb) = M_b(k,i);

        ja++; jb++; break;

      case -1:
        for(k = 0; k < M_b.extent(0); k++)
          D_b(k,jb) = M_b(k,i);

        jb++; break;

      case 0:
        break;
    }
  }
}

//Initialize the M matrix with un-backflow orbitals
void BackFlowClass::InitM(SystemClass &system)
{
  int i,j;

  for(i = 0; i < Ma.extent(0); i++)
    for (j = 0; j < Ma.extent(1); j++)
      Ma(i,j) = eigs(i,j);

  for(i = 0; i < Mb.extent(0); i++)
    for (j = 0; j < Mb.extent(1); j++)
      Mb(i,j) = eigs(i + Ma.extent(0),j);

}

void BackFlowClass::BackFlowTransform(SystemClass &system,int L)
{
  int i,k,l,site,neighbor,coeff[4];
  Ma_bf = 0; Mb_bf = 0;

#pragma omp parallel for private(i,k,l,site,neighbor,coeff)
  for(k = 0; k < Ma.extent(0); k++) // loop in state
  { 
    for(i = 0; i < BackflowNeighbors.size(); i++) // loop in site
    { 
      site     = (BackflowNeighbors[i].first).first;
      neighbor = (BackflowNeighbors[i].first).second;
      l        = BackflowNeighbors[i].second;

      coeff[0] = spin(2*site)*spin(2*site+1)*(1-spin(2*neighbor))*(1-spin(2*neighbor+1));
      coeff[1] = spin(2*site)*(1-spin(2*site+1))*spin(2*neighbor+1)*(1-spin(2*neighbor));
      coeff[2] = spin(2*site)*spin(2*site+1)*spin(2*neighbor+1)*(1-spin(2*neighbor));
      coeff[3] = spin(2*site)*(1-spin(2*site+1))*(1-spin(2*neighbor))*(1-spin(2*neighbor+1));

      Ma_bf(k,site) += ((complex<double>)coeff[0] * var(1+4*l) + (complex<double>)coeff[1] * var(2+4*l) 
                    +  ( complex<double>)coeff[2] * var(3+4*l) + (complex<double>)coeff[3] * var(4+4*l)) 
                    *  Ma(k,neighbor);
        
    }//end i
    for(i = 0; i < Ma.extent(1); i++)
      Ma_bf(k,i) += var(0) * Ma(k,i);
  } //end of the alpha electron 

#pragma omp parallel for private(i,k,l,site,neighbor,coeff)
  for(k = 0; k < Mb.extent(0); k++) // loop in state
  {
    for(i = 0; i < BackflowNeighbors.size(); i++) // loop in site
    {
      site     = (BackflowNeighbors[i].first).first;
      neighbor = (BackflowNeighbors[i].first).second;
      l        = BackflowNeighbors[i].second;

      coeff[0] = spin(2*site)*spin(2*site+1)*(1-spin(2*neighbor))*(1-spin(2*neighbor+1));
      coeff[1] = spin(2*site+1)*(1-spin(2*site))*spin(2*neighbor)*(1-spin(2*neighbor+1));
      coeff[2] = spin(2*site)*spin(2*site+1)*spin(2*neighbor)*(1-spin(2*neighbor+1));
      coeff[3] = spin(2*site+1)*(1-spin(2*site))*(1-spin(2*neighbor))*(1-spin(2*neighbor+1));

      Mb_bf(k,site) += ((complex<double>)coeff[0] * var(1+4*l) + (complex<double>)coeff[1] * var(2+4*l) 
                    +  ( complex<double>)coeff[2] * var(3+4*l) + (complex<double>)coeff[3] * var(4+4*l)) 
                    *    Mb(k,neighbor);
    }//end i
    for(i = 0; i < Mb.extent(1); i++) 
      Mb_bf(k,i) += var(0) * Mb(k,i);
  } //end the beta electron
}

void BackFlowClass::BackFlowCoord(SystemClass &system)
{
  BackFlowTransform(system,level);
}

//This is the initialize for the parameters
void BackFlowClass::GetParameters(SystemClass &system)
{
  ifstream infile;
  infile.open("parameters.txt");
 
  infile >> level;

  cout << "The iterative level: " << level << endl;

  if(level > 0)
    NumParams = 4 * level + 1;
  else
    NumParams = 0;

  if (level > 0)
    var.resize(NumParams);
 
  for(int i = 0; i < NumParams; i++)
  {
    infile >> var(i).real(); 
    infile >> var(i).imag(); 
  }

  infile.close();

  cout << "The variational parameters:" << endl;
  for(int i = 0; i < var.size(); i++)
    cout << var(i) << " "; 
  cout << endl;

  if(level == 0)
    cout << "No backflow is applied" << endl;

  if(level > 0)
  {
    infile.open("backflow_neighbors.txt");

    if (!infile){
      cout<<"Could not open backflow_neighbors.txt "<<endl;
      cout<<"Aborting!"<<endl;
      exit(1);
    }

    while (!infile.eof()){
      int i;int j;int k;
      infile>>i;infile>>j;infile>>k;
      if (!infile.eof())
        BackflowNeighbors.push_back(make_pair(make_pair(i,j),k));
    }

    infile.close();
    cout << BackflowNeighbors.size() << endl;
    cout << "End the reading of backflow_neighbors.txt" << endl;
  }
}

void BackFlowClass::PrintMatrix(Array<complex<double>,2> &M)
{
  for(int i = 0; i < M.extent(0); i++)
  {
    for(int j = 0; j < M.extent(1); j++)
      cout << M(i,j) << " ";
    cout << endl;
  }

}

void BackFlowClass::PrintMatrix(SmartEigen &Eigs, int size)
{
  for(int i = 0; i < size; i++)
  {
    for(int j = 0; j < size; j++)
      cout << Eigs.M(i,j) << " ";
    cout << endl;
  }

}

void BackFlowClass::PrintSpin(SystemClass &system)
{
  for(int i = 0; i < system.x.size(); i++)
    cout << i << " " << system.x(i) << endl;

  cout << endl;
}

complex<double> BackFlowClass::evaluate(SystemClass &system)
{
  if(level > 0)
  {
    ConvertSpin(system);
    BackFlowCoord(system);
    FillDet(system,DA,DB,Ma_bf,Mb_bf);
  }
  else
    FillDet(system,DA,DB,Ma,Mb);

  DET_a = DA.Det(); DET_b = DB.Det(); //* Sign(system);

}
complex<double> BackFlowClass::evaluateRatio(SystemClass &system,int swap1, int swap2){}
complex<double> BackFlowClass::evaluateRatio_check(SystemClass &system, int swap1, int swap2){}
void BackFlowClass::Swap(int i, int j){}
void BackFlowClass::Move(int site, int end_site, int s) {}
void BackFlowClass::UpdateDets(SystemClass &system,int swap1, int swap2){}
 
