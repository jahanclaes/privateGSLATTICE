#include "RVB_fast.h"
#include "MatrixOps.h"
// #include "FT.h"
#include "SmartMatrix.h"
#include <algorithm>

complex<double> 
RVBFastPsiClass::logevaluate(SystemClass &system,int &sign)
{
  assert(1==2);
}



double
RVBFastPsiClass::Psi_alpha_over_psi_FD(SystemClass &system)
{
}

//assumes swap1 and swap2 are of different spin
complex<double> 
RVBFastPsiClass::evaluateRatio_noStore(SystemClass &system,int swap1, int swap2)
{
  ///let's define spin up as swap1
  if (system.x(swap1)!=0){
    swap(swap1,swap2);
  }

  //evaluate the new spin up column
  
  //swap1 has been set to be the spin up value
  //loop over the spin down particles
  for (int j=0;j<system.x.size();j++)
    if (system.x(j)==1){
      u(mat.DetPos(j))=Phi(swap1,j,system);
    }

  //swap2 has been set to be the spin down value
  //loops over the spin up particles
  for (int i=0;i<system.x.size();i++)
    if (system.x(i)==0){
      up(mat.DetPos(i))=Phi(i,swap2,system);
    }

  complex<double> sp=mat.UpdateRowAndInverse(mat.DetPos(swap2),up);
  complex<double> s=mat.UpdateColAndInverse(mat.DetPos(swap1),u);
  //  if (fabs(s)<1e-10 || fabs(sp)<1e-10)
  //    mat.CalcAndSaveInverse();
  ///    cerr<<"Should rebuild"<<endl;
    
  complex<double> actRatio;
  if ((s.real()==0  && s.imag()==0)|| (sp.real()==0 && sp.imag()==0))
    actRatio=0;
  else 
    actRatio=exp(log(fabs(abs(s)))+log(fabs(abs(sp))))*(s/abs(s))*(sp/abs(sp));
  ////  cerr<<"ACTRATIO IS "<<s<<" "<<sp<<" "<<actRatio<<endl;
  //  evaluateRatio_check(system,swap1,swap2);
  /////  cerr<<"Done actratio checking"<<endl;
  return -actRatio;
}


complex<double> 
RVBFastPsiClass::evaluateRatio_check(SystemClass &system, int swap1, int swap2)
{

  SmartMatrix mat_check;
  mat_check.Init(mat.M.extent(0));
  system.Swap(swap1,swap2);
  FillDet(system,mat_check);
  complex<double> pre=mat_check.Det();
  system.Swap(swap1,swap2);
  FillDet(system,mat_check);
  complex<double> post=mat_check.Det();
  cerr<<"Check ratio: "<<post<<" "<<pre<<" "<<post/pre<<endl;
  return post/pre;
}
 





void 
RVBFastPsiClass::FillDet_check(SystemClass &system,SmartMatrix &myMat)
{
  int upDet=-1;
  int downDet=-1;
  bool ok=true;
  for (int i=0;i<system.x.size();i++){
    if (system.x(i)==0){
      upDet++;
      downDet=-1;
      for (int j=0;j<system.x.size();j++){
	if (system.x(j)==1){
	  downDet++;
	  if (!(myMat.M(myMat.DetPos(i),myMat.DetPos(j))==Phi(i,j,system))){
	    cerr<<i<<" "<<j<<" "<<myMat.DetPos(i)<<" "<<myMat.DetPos(j)<<" "<<myMat.M(myMat.DetPos(i),myMat.DetPos(j))<<" "<<Phi(i,j,system)<<" "<<endl;
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
