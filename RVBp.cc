#include "RVBp.h"
#include "MatrixOps.h"
// #include "FT.h"
#include <Eigen/Dense>
#include "SmartEigen.h"
#include <algorithm>

/////spin 0 -> first index
////spin 1 -> second index

void 
RVBpPsiClass::Init(SystemClass &system)
{
  NeedFrequentReset=true;
  //  NeedFrequentReset=false;
  cerr<<"Starting Init"<<endl;
  NumSpinUp=0;
  for (int i=0;i<system.x.size();i++){
    NumSpinUp += ((system.x(i)==1 || system.x(i)==2) ? 1: 0);
  }
  mat.Init(NumSpinUp,system.x.size());
  cerr<<"RVB BINS"<<endl;
  PairingFunction.Init(system);
  cerr<<"RVB BINS done"<<endl;
  ReadParams=false;

  if (ReadParams){
    cerr<<"READING IN THE PARMS HERE"<<endl;
    ifstream infile;
    infile.open("params.dat");
    assert(infile);
    int num=0;
    while (!infile.eof()){
      double numGarbage;
      double val;
       infile>>val;
       infile>>numGarbage;
       if (!infile.eof())
 	PairingFunction.f0[num]=val;
       num++;
     }
     infile.close();
     
     for (int i=0;i<system.rList.size();i++)
       for (int j=0;j<system.rList.size();j++){
	 int bin=PairingFunction.FindBin(i,j);
	 cerr<<"READ BIN VALUES: "<<i<<" "<<j<<" "<<PairingFunction.f0[bin].real()<<endl;
       }
   }
  
  ReadPairingFunction=false;
  if (ReadPairingFunction){
    int bin=0;
    ifstream infile;
    infile.open("PairingFunction.dat");
    assert(infile);
    while (!infile.eof()){
      int i;	int j; double tr; double ti;
	assert(!infile.eof());
	infile>>i;
	if (!infile.eof()){
	  assert(!infile.eof());
	  infile>>j;
	  bin=PairingFunction.FindBin(i,j);
	  assert(bin<PairingFunction.f0.size());
	assert(!infile.eof());
	double real; infile>>real;;
	PairingFunction.f0[bin].real(real);
	
	PairingFunction.f0[bin].real(); //=PairingFunction.f0[bin].real()/100.0;
	assert(!infile.eof());
	double imag; infile>>imag;
	PairingFunction.f0[bin].imag(imag);
	PairingFunction.f0[bin].imag(); //=PairingFunction.f0[bin].imag()/100.0;
      }
	cerr<<i<<" "<<j<<" "<<PairingFunction.f0[bin]<<" "<<bin<<endl;
    }
    infile.close();
    cerr<<"READ IN PAIRING FUNCTIONS"<<endl;
  }
  cerr<<"done writing "<<endl;
  NumParams=PairingFunction.f0.size();
  evaluate(system);
  cerr<<"Init Done"<<endl;
}



void 
RVBpPsiClass::Swap(int i, int j)
{
  swap(mat.UpPos[i],mat.UpPos[j]);
  swap(mat.DownPos[i],mat.DownPos[j]);
  //  swap(mat.DetPos[i],mat.DetPos[j]);
}


void RVBpPsiClass::Move(int site, int end_site, int spin)
{
  //  cerr<<"Moving from site "<<site<<" to "<<end_site<<" with spin "<<spin <<endl;
  if (spin==1){
    mat.UpPos[end_site]=mat.UpPos[site];
    mat.UpPos[site]=-1;
  }
  else{
    mat.DownPos[end_site]=mat.DownPos[site];
    mat.DownPos[site]=-1;
  }

}


complex<double>
RVBpPsiClass::Deriv(SystemClass &system,int bin)
{
  assert(1==2);
  complex<double> totalRatio=0.0;
  for (int i=0;i<system.x.size();i++){
    for (int j=0;j<system.x.size();j++){
      if (system.x(i)==1 && system.x(j)==-1){
	if (PairingFunction.PhiMap(i,j)==bin){
	  totalRatio+=mat.MInverse(mat.DownPos[j],mat.UpPos[i]);
	}
	else {
	  
	}
      }
    }
  }
  return totalRatio;
}


double RVBpPsiClass::GetParam_real(int i)
{
  return PairingFunction.f0[i].real();

}

double RVBpPsiClass::GetParam_imag(int i)
{
  return PairingFunction.f0[i].imag();

}


void RVBpPsiClass::SetParam_real(int i, double param)
{
  PairingFunction.f0[i].real(param);

}
void RVBpPsiClass::SetParam_imag(int i, double param)
{
  PairingFunction.f0[i].imag(param);

}



void
RVBpPsiClass::AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs)
{
  AllDerivs(system, derivs,0,derivs.size());

}

void
RVBpPsiClass::AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start,int stop)
{

  //  mat.CheckInverse();
  
  for (int i=start;i<stop;i++)
    derivs(i)=0.0;
  for (int i=0;i<system.x.size();i++){
    for (int j=0;j<system.x.size();j++){
      if ((system.x(i)==1 || system.x(i)==2)  && (system.x(j)==-1 || system.x(j)==2)){
	int bin=start+PairingFunction.PhiMap(i,j);
	//	cerr<<"I think the deriv for "<<bin<<" "<<i<<" "<<j<<" "<<system.x(i)<<" "<<system.x(j)<<" "<<" is "<<mat.MInverse(mat.DownPos[j],mat.UpPos[i])<<" "<<TestDerivs(PairingFunction.PhiMap(i,j),system)<<endl;
	if (system.x(i)==2 && system.x(j)==2)
	  derivs(bin)+=mat.MInverse(mat.DownPos[j],mat.UpPos[i]); // +mat.MInverse(mat.DownPos[i],mat.UpPos[j]);
	else
	  derivs(bin)+=mat.MInverse(mat.DownPos[j],mat.UpPos[i]);
      }
      else {
	
      }
    }
  }

//   for (int i=0;i<system.x.size();i++){
//     for (int j=0;j<system.x.size();j++){
//       if ((system.x(i)==1 || system.x(i)==2)  && (system.x(j)==-1 || system.x(j)==2)){
// 	int bin=start+PairingFunction.PhiMap(i,j);
// 	if (abs(derivs(bin)-TestDerivs(PairingFunction.PhiMap(i,j),system))>1e-5){
	  
// 	  cerr<<"I think the deriv for "<<bin<<" "<<i<<" "<<j<<" "<<system.x(i)<<" "<<system.x(j)<<" "<<" is "<<derivs(bin)<<" "<<TestDerivs(PairingFunction.PhiMap(i,j),system)<<" "<<abs(derivs(bin)-TestDerivs(PairingFunction.PhiMap(i,j),system))<<endl;
// 	}
//       }
//     }
//   }
}
double RVBpPsiClass::TestDerivs(int derivInt,SystemClass &system)
{
  //  cerr<<"Testing "<<derivInt<<endl;
  double currParam=GetParam_real(derivInt);
//   for (double myStep=-0.01;myStep<0.01;myStep+=0.001){
//     double step_energy=0.0;
//     double countSteps=0.0;
//     SetParam_real(derivInt,currParam+myStep);
//     cerr<<myStep<<" "<<evaluate_noInverse(system)<<endl;
//   }
  
  SetParam_real(derivInt,currParam+0.001);
  double up=evaluate_noChange(system).real();
  SetParam_real(derivInt,currParam-0.001);
  double down=evaluate_noChange(system).real();

  SetParam_real(derivInt,currParam);
  double curr=evaluate_noChange(system).real();
  //  cerr<<"Done Testing "<<derivInt<<endl;
  return  (up-down)/(0.002*curr);
}


complex<double>
RVBpPsiClass::Phi(int i,int j,SystemClass &system)
{
  return PairingFunction.Phi(i,j);
  double gutz_correct=PairingFunction.gutz(i,j,system);
  if (fabs(gutz_correct-PairingFunction.Phi(i,j).real())>1e-10){
    cerr<<"ERROR: "<<gutz_correct<<" "<<PairingFunction.Phi(i,j)<<" "<<i<<" "<<j<<endl;
  }
  complex<double> toReturn(gutz_correct,0.0);
  return toReturn;
}

void RVBpPsiClass::SetParams(double delta3,SystemClass &system)

{
  SetParams(0,delta3,system);
}

void 
RVBpPsiClass::SetParams(int i,double delta3, SystemClass &system)
{
  assert(1==2);
}


complex<double> 
RVBpPsiClass::evaluate_noChange(SystemClass &system)
{
  SmartEigen mat_check;
  mat_check.Init(mat.M.rows(),mat.UpPos.size());
  FillDet(system,mat_check);
  return mat_check.Det();
}

complex<double> 
RVBpPsiClass::evaluate_noInverse(SystemClass &system)
{
  FillDet(system,mat);
  return mat.Det();
}


complex<double>
RVBpPsiClass::evaluate(SystemClass &system)
{
  complex<double> myAns=evaluate_noInverse(system);
  mat.CalcAndSaveInverse();
  return myAns;
}


void
RVBpPsiClass::Reject(SystemClass &system,int swap1,int swap2)
{
  
}


void RVBpPsiClass::Reject(SystemClass &system,int site,int end_site,int spin)
{


}

void 
RVBpPsiClass::UpdateDets(SystemClass &system,int swap1, int swap2)
{
  mat.InverseUpdate(colIndices,rowIndices,newColsp,newRowsp);  
}

void 
RVBpPsiClass::UpdateDets(SystemClass &system,int site, int end_site,int spin)
{

  if (spin==1){
    mat.UpdateRowAndInverse(mat.UpPos[end_site],col);
  }
  else{
    mat.UpdateColAndInverse(mat.DownPos[end_site],col);
  }


}

complex<double>
RVBpPsiClass::evaluateRatio(SystemClass &system,int start, int stop, int spin)
{

  int countParity=0;
  int myMin=min(start,stop);
  int myMax=max(start,stop);
  if (spin==1){
    for (int i=myMin+1;i<myMax;i++)
      if (mat.UpPos[i]!=-1)
	countParity++;
  }
  else {
    for (int i=myMin+1;i<myMax;i++)
      if (mat.DownPos[i]!=-1)
	countParity++;
  }
  countParity = ( (countParity % 2) ==0 ) ?  1: -1;
  int not_spin = (spin==1) ? -1 : 1;
  col=Eigen::VectorXcd::Zero(NumSpinUp);
  for (int site=0;site<system.x.size();site++){
    //    cerr<<"Site vals: "<<system.x(site)<<endl;
    if (spin==1 && (system.x(site)==not_spin ||  system.x(site)==2)){
      //      cerr<<"UP: "<<col.size()<<" "<<mat.DownPos[site]<<" "<<mat.M.rows()<<" "<<mat.M.cols()<<endl;
      col(mat.DownPos[site])=Phi(stop,site,system);
    }
    else if (spin==-1 && (system.x(site)==not_spin ||  system.x(site)==2)){
      col(mat.UpPos[site])=Phi(site,stop,system);
    }
  }
  if (spin==1){
    complex<double> ratio= mat.RowRatio(mat.UpPos[stop],col);
    //    complex<double> ratio_check=evaluateRatio_check(system,start,stop,spin);
    //    cerr<<"GRR: "<<ratio<<" "<<ratio_check<<endl;
    ratio.real(ratio.real()*countParity);
    ratio.imag(ratio.imag()*countParity);
    //    ratio.real()*=countParity;
    //    ratio.imag()*=countParity;
    return ratio; //_check;
  }
  else{

    complex<double> ratio= mat.ColRatio(mat.DownPos[stop],col);
    //    complex<double> ratio_check=evaluateRatio_check(system,start,stop,spin);
    //    cerr<<"GRR: "<<ratio<<" "<<ratio_check<<endl;
    ratio.real(ratio.real()*countParity);
    ratio.imag(ratio.imag()*countParity);
    //    ratio.real()*=countParity;
    //    ratio.imag()*=countParity;
    //    cerr<<"Ratios B are: "<<ratio<<" "<<ratio_check<<endl;

    return ratio;// _check;
  }
}

complex<double> 
RVBpPsiClass::evaluateRatio_check(SystemClass &system, int site, int end_site,
				     int spin)
{

  SmartEigen mat_check;
  mat_check.Init(mat.M.rows(),mat.UpPos.size());
  system.Move(end_site,site,spin);
  FillDet(system,mat_check);
  complex<double> pre=mat_check.Det();
  system.Move(site,end_site,spin);
  FillDet(system,mat_check);
  //  cerr<<mat_check.M<<endl;
  //  cerr<<endl;
  //  cerr<<endl;
  //  cerr<<mat.M-mat_check.M<<endl;
  complex<double> post=mat_check.Det();
  return post/pre;
}


//assumes swap1 and swap2 are of different spin
complex<double> 
RVBpPsiClass::evaluateRatio(SystemClass &system,int swap1, int swap2)
{

  int maxSwap=max(swap1,swap2);
  int minSwap=min(swap1,swap2);
  int mySign=1;
  if ((maxSwap-minSwap) % 2==0)
    mySign=-1;
  else
    mySign=1;

  //  mat.SaveInverse();
  ///let's define spin up as swap1
  if (system.x(swap1)!=1){
    swap(swap1,swap2);
  }

  //evaluate the new spin up column
  newRowsp.resize(1,NumSpinUp);
  newColsp.resize(NumSpinUp,1);

  colIndices.resize(1);
  colIndices[0]=mat.UpPos[swap1];
  rowIndices.resize(1);
  rowIndices[0]=mat.DownPos[swap2];

  //swap1 has been set to be the spin up value
  //loop over the spin down particles
  for (int j=0;j<system.x.size();j++){
    if ((system.x(j)==-1) || ((system.x(j)==2))){  
      newColsp(mat.DownPos[j],0)=Phi(swap1,j,system);
    }
    //swap2 has been set to be the spin down value
    //loops over the spin up particles
    if ((system.x(j)==1) || (system.x(j)==2)){
      newRowsp(0,mat.UpPos[j])=Phi(j,swap2,system);
    }
  }

  complex<double> test_ratio=mat.Ratio_ncol_nrowp(colIndices,rowIndices,newColsp,newRowsp);

  int countParity=0;
  int myMin=min(swap1,swap2);
  int myMax=max(swap1,swap2);
  for (int i=myMin+1;i<myMax;i++){
    //    if ( (mat.UpPos[i]!=-1)  || (mat.DownPos[i]!=-1))//HACK!
    if ( (mat.UpPos[i]!=-1) )
      countParity++;
    if (mat.DownPos[i]!=-1)
      countParity++;
  }
  test_ratio*= (( (countParity % 2)==0)  ? 1 :-1);

  //  complex<double> check_ratio = evaluateRatio_check(system,swap1,swap2);
  //  complex<double> diff=test_ratio-check_ratio;
  //  cerr<<"CHECKING: "<<diff<<" "<<check_ratio<<" "<<test_ratio<<endl;
  //  assert(fabs(diff)<1e-4);
  rebuild=false;
  //  cerr<<"The ratio is "<<test_ratio<<endl;
  
  return test_ratio;
}



void 
RVBpPsiClass::FillDet(SystemClass &system,SmartEigen &myMat)
{
  //  cerr<<"FILLING NOW "<<system.x.size()<<endl;
  //  for (int i=0;i<system.x.size();i++){
  //    cerr<<system.x(i)<<endl;
  //  }
  for (int i=0;i<myMat.UpPos.size();i++){
    myMat.UpPos[i]=-1;
    myMat.DownPos[i]=-1;
  }
  int upDet=-1;
  int downDet=-1;
  for (int i=0;i<system.x.size();i++){
    if ( (system.x(i)==1) || (system.x(i)==2) ){ //HACK!
      upDet++;
      downDet=-1;
      for (int j=0;j<system.x.size();j++){
	if ((system.x(j)==-1) || (system.x(j)==2) ){ //HACK!
	  downDet++;
	  //	  cerr<<"GRR GRR"<<upDet<<" "<<downDet<<endl;
	  //	  cerr<<"The phi is "<<i<<" "<<j<<" "<<Phi(i,j,system)<<endl;
	  myMat.M(upDet,downDet)=Phi(i,j,system);
	  myMat.UpPos[i]=upDet;
	  myMat.DownPos[j]=downDet;
	  //	  myMat.DetPos[i]=upDet;
	  //	  myMat.DetPos[j]=downDet;
	}
      }
    }
  }
  ///  cerr<<"WRITING ME "<<endl;
  //  for (int i=0;i<myMat.M.rows();i++){
  //    for (int j=0;j<myMat.M.cols();j++){
  //      cerr<<myMat.M(i,j)<<" ";
  //    }
  //    cerr<<endl;
  //  }
  //  cerr<<"DONE WRITING ME "<<endl;
}

complex<double> 
RVBpPsiClass::evaluateRatio_check(SystemClass &system, int swap1, int swap2)
{

  SmartEigen mat_check;
  mat_check.Init(mat.M.rows());
  system.Swap(swap1,swap2);
  FillDet(system,mat_check);
  complex<double> pre=mat_check.Det();
  system.Swap(swap1,swap2);
  FillDet(system,mat_check);
  complex<double> post=mat_check.Det();
  //  cerr<<"Check ratio: "<<post<<" "<<pre<<" "<<post/pre<<endl;
  return post/pre;
}

