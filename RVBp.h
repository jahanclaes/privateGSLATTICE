///TEST TEST TEST
#ifndef RVBP_H
#define RVBP_H
#include <iostream>
#include "Blitz.h"
#include "Random/Random.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include "SystemClass.h"
#include "WaveFunction.h"
#include "SmartEigen.h"
#include "SmartMatrix.h"
#include "PairingFunctionAllBin.h"


class RVBpPsiClass  : public WaveFunctionClass
{
public:
  SmartEigen mat;
  SmartMatrix matp;
  bool ReadParams;
  bool ReadPairingFunction;
  //  RVBpPsiClass::
 RVBpPsiClass(PairingFunctionAllBin &pf)  : PairingFunction(pf) {
    ReadParams=false;
    ReadPairingFunction=false;
    Name="RVBpPsi";
  }

  double TestDerivs(int derivInt,SystemClass &system);
  complex<double> evaluate_noChange(SystemClass &system);
  
  complex<double> evaluateRatio(SystemClass &system,int start, int stop, int spin);
  vector<Array<complex<double>  ,1> > newCols;
  vector<Array<complex<double>  ,1> > newRows;
  Eigen::MatrixXcd newColsp;
  Eigen::MatrixXcd newRowsp;

  vector<int> colIndices;
  vector<int> rowIndices;
  void Copy(WaveFunctionClass* wf)
  {
    RVBpPsiClass &b(*(RVBpPsiClass*)wf);
    mat=b.mat;
    newCols=b.newCols;
    newRows=b.newRows;
    colIndices=b.colIndices;
    rowIndices=b.rowIndices;
    NumSpinUp=b.NumSpinUp;
    rebuild=b.rebuild;
    

  }

  PairingFunctionAllBin &PairingFunction; 

  int NumSpinUp;
  bool rebuild;
  Array<complex<double>,1> u;  Array<complex<double>,1> up;
  void SetParams(int i,double delta3, SystemClass &system);

  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start,int stop);
  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs);

  double GetParam_real(int i);
  void SetParam_real(int i, double param);

  double GetParam_imag(int i);
  void SetParam_imag(int i, double param);

  complex<double> Deriv(SystemClass &system,int bin); 

  void Init(SystemClass &system);
  void Swap(int i, int j);
  void Move(int site, int end_site,int spin);
  void SetParams(double delta3,SystemClass &system);
  
  complex<double>Phi(int i,int j,SystemClass &system);

  void FillDet(SystemClass &system, SmartEigen &myMat);
  void FillDet(SystemClass &system, SmartMatrix &myMat);
  complex<double> evaluate(SystemClass &system);
  complex<double> evaluate_noInverse(SystemClass &system);
  complex<double> evaluateRatio(SystemClass &system,int swap1, int swap2);
  complex<double> evaluateRatio_energy(SystemClass &system,int swap1, int swap2);
  complex<double>  evaluateRatio_check(SystemClass &system, int site, int end_site,
				       int spin);

  double Sign(SystemClass &system);

  void Reject(SystemClass &system,int site,int end_site,int spin);
  complex<double> evaluateRatio_swap(SystemClass &system, int swap1, int swap2);

  void Reject(SystemClass &system,int swap1,int swap2);
  void UpdateDets(SystemClass &system,int swap1, int swap2);
  void UpdateDets(SystemClass &system,int site, int end_site,int spin);

  Eigen::VectorXcd col; 
  Array<complex<double> ,1> colp; 
  complex<double> logevaluate(SystemClass &system,int &sign)
    {
      assert(1==2);
    }
  
  std::complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2);
  




void  FillDet_check(SystemClass &system,SmartEigen &myMat)
{
  int upDet=-1;
  int downDet=-1;
  bool ok=true;
  for (int i=0;i<system.x.size();i++){
    if ((system.x(i)==1) || (system.x(i)==2) ){
      upDet++;
      downDet=-1;
      for (int j=0;j<system.x.size();j++){
	if ( (system.x(j)==-1) || (system.x(j)==2)){
	  downDet++;
	  if (!(myMat.M(myMat.UpPos[i],myMat.DownPos[j])==Phi(i,j,system))){
	    cerr<<i<<" "<<j<<" "<<myMat.UpPos[i]<<" "<<myMat.DownPos[j]<<" "<<myMat.M(myMat.UpPos[i],myMat.DownPos[j])<<" "<<Phi(i,j,system)<<" "<<endl;
	    ok=false;
	  }

	}

      }
    }
  }
  for (int i=0;i<system.x.size();i++)
    cerr<<system.x(i)<<" ";
  cerr<<endl;
  assert(ok);
}


/*   //For debugging */
/*   SmartEigen mat_check; */
/*   complex<double> logevaluate(SystemClass &system,int &sign); */
/*   complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2); */

/*   double Psi_alpha_over_psi_FD(SystemClass &system); */
/*   complex<double> evaluateRatio_noStore(SystemClass &system,int swap1, int swap2); */
  


/*   void FillDet_check(SystemClass &system,SmartEigen &myMat); */

};

#endif
