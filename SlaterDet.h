///TEST TEST TEST
#ifndef SLATERDET_H
#define SLATERDET_H
#include <iostream>
#include "Blitz.h"
#include "Random/Random.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include "SystemClass.h"
#include "WaveFunction.h"
#include "SmartEigen.h"
#include "PairingFunctionAllBin.h"
#include "SharedEigs.h"

class SlaterDetPsiClass  : public WaveFunctionClass
{
public:
  SmartEigen mat;
  bool ReadParams;
  bool ReadPairingFunction;

  Eigen::VectorXcd col;   


  vector<Array<complex<double>  ,1> > newCols;
  vector<Array<complex<double>  ,1> > newRows;
  vector<int> colIndices;
  vector<int> rowIndices;
  int NumSpinUp;
  bool rebuild;
  int mySpin;
  //  int NumElectrons;



  SharedEigsClass &SharedEigs;
  //  PairingFunctionAllBin &PairingFunction; 

  void RebuildParams()
  {

    //    Eigen::ColPivHouseholderQR<Eigen::MatrixXcd>  QR(SharedEigs.eigs.transpose());
    Eigen::HouseholderQR<Eigen::MatrixXcd>  QR(SharedEigs.eigs.transpose());
    //k    SharedEigs.eigs=QR.matrixQ().transpose();
    SharedEigs.eigs=QR.householderQ().transpose();
  } 

  void Copy(WaveFunctionClass* wf)
  {
    SlaterDetPsiClass &b(*(SlaterDetPsiClass*)wf);
    mat=b.mat; 
    newCols=b.newCols;
    newRows=b.newRows; 
    colIndices=b.colIndices; 
    rowIndices=b.rowIndices;
    NumSpinUp=b.NumSpinUp; 
    rebuild=b.rebuild; 
    mySpin=b.mySpin;
  }

 SlaterDetPsiClass(SharedEigsClass &e)  : SharedEigs(e) {
    ReadParams=false;
    ReadPairingFunction=false;
    Name="SlaterDet";
  }

  double TestDerivs(int derivInt,SystemClass &system);

  complex<double> evaluate_noChange(SystemClass &system);
  complex<double> evaluateRatio(SystemClass &system,int start, int stop, int spin);





  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start,int stop);
  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs);

  double GetParam_real(int i);
  void SetParam_real(int i, double param);
  double GetParam_imag(int i);
  void SetParam_imag(int i, double param);

  complex<double> Deriv(SystemClass &system,int bin); 

  void Init(SystemClass &system,int mySpin);
  void Swap(int i, int j);
  void Move(int site, int end_site,int spin);
  //  void SetParams(double delta3,SystemClass &system);
  //  void SetParams(int i,double delta3, SystemClass &system);  
  complex<double>Phi(int i,int j,SystemClass &system);

  void FillDet(SystemClass &system, SmartEigen &myMat);
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


  complex<double> logevaluate(SystemClass &system,int &sign)
    {
      assert(1==2);
    }
  
  std::complex<double> evaluateRatio_check(SystemClass &system, int a, int b)
    {
      assert(1==2);
    }


  void  FillDet_check(SystemClass &system,SmartEigen &myMat);


/*   //For debugging */
/*   SmartEigen mat_check; */
/*   complex<double> logevaluate(SystemClass &system,int &sign); */
/*   complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2); */

/*   double Psi_alpha_over_psi_FD(SystemClass &system); */
/*   complex<double> evaluateRatio_noStore(SystemClass &system,int swap1, int swap2); */
  


/*   void FillDet_check(SystemClass &system,SmartEigen &myMat); */

};

#endif
