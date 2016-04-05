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
    //  NumSpinUp=system.x.size()/2;
  mat.Init(NumSpinUp,system.x.size());
  //  u.resize(NumSpinUp);  
  //  up.resize(NumSpinUp);
  cerr<<"RVB BINS"<<endl;
  PairingFunction.Init(system);
  cerr<<"RVB BINS done"<<endl;
  ReadParams=false;

//   bool param2=false;
//   if (ReadParams){
//       cerr<<"Now reading parameters"<<endl;
//       ifstream infile;

//       infile.open("paramsAA.txt");
// 	//HACK	infile.open("paramspAA.txt");
//      dVec zero(0.0,0.0);
//      while (!infile.eof()){
//        dVec r;
//        infile>>r[0];
//        infile>>r[1];
//        //       double minDist=system.minDist(r,zero);

//        dVec minDiff=system.minDiff(r,zero);
//        int theBin=PairingFunctionAA.FindBin(minDiff);
// 	if (theBin==-1){
// 	bool found=false;
// 	for (int ip=-1;ip<=1;ip++)
// 	  for (int jp=-1;jp<=1;jp++){
// 	    dVec diffp(minDiff[0]+system.a1[0]*system.N1*ip+system.a2[0]*jp*system.N2,minDiff[1]+system.a1[1]*ip*system.N1+system.a2[1]*jp*system.N2);
// 	    dVec negDiff=-diffp;
// 	    if (theBin==-1)
// 	      theBin=PairingFunctionAA.FindBin(diffp);
// 	    if (theBin==-1)
// 	      theBin=PairingFunctionAA.FindBin(negDiff);
// 	    //	    cerr<<"Bin not found but now: "<<bin<<" "<<diff<<endl;
// 	  }
// 	}
//        if (theBin==-1){
// 	 cerr<<"Min dist is "<<minDiff<<endl;
// 	 for (int i=0;i<PairingFunctionAA.binVals.size();i++)
// 	   cerr<<"Bin Vals "<<i<<" "<<PairingFunctionAA.binVals[i]<<endl;

//        }
//        if (theBin!=-1){
// 	 cerr<<"Going to be Setting AA "<<theBin<<" to "<<PairingFunctionAA.f0[theBin]<<endl;
// 	 infile>>PairingFunctionAA.f0[theBin].real();
// 	 infile>>PairingFunctionAA.f0[theBin].imag();
// 	 cerr<<"Setting AA "<<theBin<<" to "<<PairingFunctionAA.f0[theBin]<<endl;
//        }
//      }
//      cerr<<"done with AA"<<endl;
//      infile.close();
//      if (!param2)
//        infile.open("paramsBB.txt");
//      else 
//        infile.open("paramsBB.txt");
//        //HACK       infile.open("paramspBB.txt");
//      while (!infile.eof()){
//        dVec r;
//        infile>>r[0];
//        infile>>r[1];
//        //       double minDist=system.minDist(r,zero);
//        dVec minDiff=system.minDiff(r,zero);
//        int theBin=PairingFunctionBB.FindBin(minDiff);
// 	if (theBin==-1){
// 	bool found=false;
// 	for (int ip=-1;ip<=1;ip++)
// 	  for (int jp=-1;jp<=1;jp++){
// 	    dVec diffp(minDiff[0]+system.a1[0]*system.N1*ip+system.a2[0]*jp*system.N2,minDiff[1]+system.a1[1]*ip*system.N1+system.a2[1]*jp*system.N2);
// 	    dVec negDiff=-diffp;
// 	    if (theBin==-1)
// 	      theBin=PairingFunctionBB.FindBin(diffp);
// 	    if (theBin==-1)
// 	      theBin=PairingFunctionBB.FindBin(negDiff);
// 	    //	    cerr<<"Bin not found but now: "<<bin<<" "<<diff<<endl;
// 	  }
// 	}

    	


//  	if (theBin!=-1){
//  	  infile>>PairingFunctionBB.f0[theBin].real();
//  	  infile>>PairingFunctionBB.f0[theBin].imag();
//  	  cerr<<"Min dist for BB is "<<minDiff<<endl;
//  	  cerr<<"Setting BB "<<theBin<<" to "<<PairingFunctionBB.f0[theBin]<<endl;
//  	}
//       }
//       infile.close();
//       if (!param2)
//         infile.open("paramsAB.txt");
//       else 
//         //HACK       infile.open("paramspAB.txt");
//         infile.open("paramsAB.txt");
//       while (!infile.eof()){
//         dVec r;
//         infile>>r[0];
//         infile>>r[1];
//         dVec minDiff=system.minDiff(r,zero);
//         int theBin=PairingFunctionAB.FindBin(minDiff);

//  	if (theBin==-1){
//  	bool found=false;
//  	for (int ip=-1;ip<=1;ip++)
//  	  for (int jp=-1;jp<=1;jp++){
//  	    dVec diffp(minDiff[0]+system.a1[0]*system.N1*ip+system.a2[0]*jp*system.N2,minDiff[1]+system.a1[1]*ip*system.N1+system.a2[1]*jp*system.N2);
//  	    dVec negDiff=-1.0*diffp;
//  	    if (theBin==-1)
//  	      theBin=PairingFunctionAB.FindBin(diffp);
//  	    if (theBin==-1)
//  	      theBin=PairingFunctionAB.FindBin(negDiff);
//  	    //	    cerr<<"Bin not found but now: "<<bin<<" "<<diff<<endl;
//  	  }
//  	}

//  	if (theBin!=-1){
	  
//  	  infile>>PairingFunctionAB.f0[theBin].real();
//  	  infile>>PairingFunctionAB.f0[theBin].imag();
//  	  cerr<<"Setting AB "<<theBin<<" to "<<PairingFunctionAB.f0[theBin]<<endl;
     
//         if (theBin==-1){
//  	 cerr<<"Min dist is "<<minDiff<<endl;
//  	 for (int i=0;i<PairingFunctionAB.binVals.size();i++)
//  	   cerr<<"Bin Vals "<<i<<" "<<PairingFunctionAB.binVals[i]<<endl;

//         }
//  	}


//       }
//       infile.close();



//      }



   if (ReadParams){
     cerr<<"READING IN THE PARMS HERE"<<endl;
     //    exit(1);
     ifstream infile;
     infile.open("params.dat");
     assert(infile);
     int num=0;
     while (!infile.eof()){
       //       int num;
       double numGarbage;
       double val;
       string garbage;
       //       infile>>garbage;
       //       infile>>num;
       infile>>val;
       //       infile>>numGarbage;
       infile>>numGarbage;
       //HACK!      cerr<<"Setting num val "<<num<<" "<<val<<endl;
       if (!infile.eof())
 	PairingFunction.f0[num]=val;
       num++;
 	//oldREAD	if (!infile.eof() && num>=2352)
 	//oldREAD	  PairingFunction.f0[num-2352]=val;
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
	infile>>PairingFunction.f0[bin].real();
	PairingFunction.f0[bin].real(); //=PairingFunction.f0[bin].real()/100.0;
	assert(!infile.eof());
	infile>>PairingFunction.f0[bin].imag();
	PairingFunction.f0[bin].imag(); //=PairingFunction.f0[bin].imag()/100.0;
      }
      cerr<<i<<" "<<j<<" "<<PairingFunction.f0[bin]<<" "<<bin<<endl;
    }
    infile.close();
    cerr<<"READ IN PAIRING FUNCTIONS"<<endl;
  }

//   for (int i=0;i<system.x.size();i++)
//     for (int j=0;j<system.x.size();j++){
//       int bin=PairingFunction.FindBin(i,j);
//       cerr<<"GRR: "<<i<<" "<<j<<" "<<PairingFunction.f0[bin].real()<<" "<<PairingFunction.f0[bin].imag()<<endl;
// 	//<<bin<<endl;
//     }
  cerr<<"done writing "<<endl;
  //  exit(1);
  NumParams=PairingFunction.f0.size();
  evaluate(system);
  cerr<<"Init Done"<<endl;
}


//currently only the size matters
void 
RVBpPsiClass::Swap(int i, int j)
{
  swap(mat.UpPos[i],mat.UpPos[j]);
  swap(mat.DownPos[i],mat.DownPos[j]);

  swap(mat.DetPos[i],mat.DetPos[j]);
  cerr<<"UP :";
  for (int ii=0;ii<mat.UpPos.size();ii++){
    cerr<<mat.UpPos[ii]<<" ";
  }
  cerr<<endl;

  cerr<<"Down :";
  for (int ii=0;ii<mat.DownPos.size();ii++){
    cerr<<mat.DownPos[ii]<<" ";
  }
  cerr<<endl;

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
  PairingFunction.f0[i].real()=param;

}
void RVBpPsiClass::SetParam_imag(int i, double param)
{
  PairingFunction.f0[i].imag()=param;

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
  //  ofstream infile;

//   infile.open("PF.dat");
//   int size=system.x.size();
//   for (int a=0;a<size;a++){
//     for (int b=0;b<size;b++){
//       infile<<PairingFunction.Phi(a,b).real()<<" ";
//     }
//     infile<<endl;
//   }
//   infile.close();
  // exit(1);
  //  cerr<<"i j "<<i<<" "<<j<<endl;
  //  cerr<<"The pairing function is "<<PairingFunction.Phi(i,j)<<endl;
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
  //  evaluate(system);
}

void 
RVBpPsiClass::UpdateDets(SystemClass &system,int site, int end_site,int spin)
{

  if (spin==1){
    //    cerr<<"I'm updating "<<mat.UpPos[end_site]<<endl;
    mat.UpdateRowAndInverse(mat.UpPos[end_site],col);
    //    cerr<<"I'm done updating"<<endl;
  }
  else{
    //    cerr<<"I'm d updating "<<mat.DownPos[end_site]<<endl;
    mat.UpdateColAndInverse(mat.DownPos[end_site],col);
    //    cerr<<"I'm done updating"<<endl;
  }


}

complex<double>
RVBpPsiClass::evaluateRatio(SystemClass &system,int start, int stop, int spin)
{
//   cerr<<"In evaluate ratio "<<endl;
//   for (int i=0;i<mat.UpPos.size();i++){
//     cerr<<i<<" "<<mat.UpPos[i]<<endl;
//   }
//   cerr<<endl;
//   for (int i=0;i<mat.DownPos.size();i++){
//     cerr<<i<<" "<<mat.DownPos[i]<<endl;
//   }

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
  //  cerr<<"DIFF IS "<<countParity<<" "<<system.CountElectrons(start,stop,spin)<<endl;
  countParity = ( (countParity % 2) ==0 ) ?  1: -1;
  //
  int not_spin = (spin==1) ? -1 : 1;
  col=Eigen::VectorXcd::Zero(NumSpinUp);
  for (int site=0;site<system.x.size();site++){
    //    cerr<<"Site vals: "<<system.x(site)<<endl;
    if (spin==1 && (system.x(site)==not_spin ||  system.x(site)==2)){
      //      cerr<<"UP: "<<col.size()<<" "<<mat.DownPos[site]<<" "<<mat.M.rows()<<" "<<mat.M.cols()<<endl;
      col(mat.DownPos[site])=Phi(stop,site,system);
    }
    else if (spin==-1 && (system.x(site)==not_spin ||  system.x(site)==2)){
      //      cerr<<"DOWN: "<<col.size()<<" "<<mat.UpPos[site]<<endl;
      col(mat.UpPos[site])=Phi(site,stop,system);
    }
  }
  //  cerr<<"HERE"<<endl;
  if (spin==1){
    complex<double> ratio= mat.RowRatio(mat.UpPos[stop],col);
    //    cerr<<"A"<<endl;
    //    complex<double> ratio_check=evaluateRatio_check(system,start,stop,spin);
    //    cerr<<"B"<<endl;
    ratio.real()*=countParity;
    ratio.imag()*=countParity;
    //    cerr<<"Ratios A are: "<<ratio<<" "<<ratio_check<<endl;

    return ratio; //_check;
  }
  else{

    complex<double> ratio= mat.ColRatio(mat.DownPos[stop],col);
    //    cerr<<"AA"<<endl;

    //    complex<double> ratio_check=evaluateRatio_check(system,start,stop,spin);
    //    cerr<<"BB"<<endl;
    ratio.real()*=countParity;
    ratio.imag()*=countParity;

    //    cerr<<"Ratios B are: "<<ratio<<" "<<ratio_check<<endl;
    return ratio;// _check;
  }
}

complex<double> 
RVBpPsiClass::evaluateRatio_check(SystemClass &system, int site, int end_site,
				     int spin)
{
  //  cerr<<"updating "<<endl;
  //  UpdateDets(system,site,end_site,spin);
  //  cerr<<"checking"<<endl;
  //FillDet_check(system,mat);
  //  assert(1==2);

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
  if (system.x(swap1)!=-1){
    swap(swap1,swap2);
  }

  //evaluate the new spin up column
  newCols.resize(1);
  newRows.resize(1);
  //  newCols[0].resize(u.size());
  //  newRows[0].resize(u.size());
  newCols[0].resize(NumSpinUp);
  newRows[0].resize(NumSpinUp);
  
  newRowsp.resize(1,NumSpinUp);
  newColsp.resize(NumSpinUp,1);

  colIndices.resize(1);
  //  colIndices[0]=mat.DownPos[swap1];
  colIndices[0]=mat.DetPos[swap1];
  assert(mat.DetPos[swap1]==mat.UpPos[swap1]);
  //  cerr<<"UP: "<<mat.DownPos[swap1]<<endl;

  rowIndices.resize(1);
  //  rowIndices[0]=mat.UpPos[swap2];
  rowIndices[0]=mat.DetPos[swap2];
  assert(mat.DetPos[swap2]==mat.DownPos[swap2]);
  //  cerr<<"UP: "<<mat.UpPos[swap2]<<endl;
  //  cerr<<"UP: "<<mat.DownPos[swap2]<<endl;

  //swap1 has been set to be the spin up value
  //loop over the spin down particles
  for (int j=0;j<system.x.size();j++){
    if (system.x(j)==1){ 
      //      u(mat.DetPos(j))=Phi(swap1,j,system);
      //      newCols[0](mat.DetPos[j])=Phi(swap1,j,system);
      //      cerr<<"A "<<j<<" "<<mat.DownPos[j]<<endl;
      //      newColsp(mat.DownPos[j],0)=Phi(swap1,j,system);
      newColsp(mat.DetPos[j],0)=Phi(swap1,j,system);
    }

  //swap2 has been set to be the spin down value
  //loops over the spin up particles
  //  for (int i=0;i<system.x.size();i++)
    else if (system.x(j)==-1){
      //      up(mat.DetPos(i))=Phi(i,swap2,system);
      //      newRows[0](mat.DetPos[j])=Phi(j,swap2,system);
      //      cerr<<"B "<<j<<" "<<mat.UpPos[j]<<endl;
      //      newRowsp(0,mat.UpPos[j])=Phi(j,swap2,system);
      newRowsp(0,mat.DetPos[j])=Phi(j,swap2,system);
    }
  }

  complex<double> test_ratio=mat.Ratio_ncol_nrowp(colIndices,rowIndices,newColsp,newRowsp);
  complex<double> check_ratio = evaluateRatio_check(system,swap1,swap2);
  int countParity=0;
  int myMin=min(swap1,swap2);
  int myMax=max(swap1,swap2);
  for (int i=myMin+1;i<myMax;i++)
    if ( (mat.UpPos[i]!=-1)  || (mat.DownPos[i]!=-1))
      countParity++;


  //  double parity=system.CountElectrons(minSwap,maxSwap,1)+system.CountElectrons(minSwap,maxSwap,-1);
  test_ratio*= (( (countParity % 2)==0)  ? 1 :-1);
  complex<double> diff=test_ratio-check_ratio;
  cerr<<"CHECKING: "<<diff<<" "<<check_ratio<<" "<<test_ratio<<endl;
      //HACK!      assert((diff*conj(diff)).real()<1e-10);
  //  cerr<<"test ratio is "<<test_ratio<<" "<<endl; //evaluateRatio_check(system,swap1,swap2)<<endl;

  rebuild=false;
  //  cerr<<"CURRENT RATIO IS "<<test_ratio<<endl;
  //  test_ratio.real()=-1*test_ratio.real(); //*mySign;
  //  test_ratio.imag()=-1*test_ratio.imag(); //*mySign;
  return test_ratio;
  //return check_ratio;
}

double 
RVBpPsiClass::Sign(SystemClass &system)
{
  assert(1==2);
    Array<int,1> detPos(mat.DetPos.size());
    for (int i=0;i<mat.DetPos.size();i++)
      if (system.x(i)==0)
	detPos(i)=mat.DetPos[i];
      else 
	detPos(i)=mat.DetPos[i]+system.x.size()/2;
    int numSwap=0;
    for (int i=0;i<detPos.size();i++){
      while (detPos(i)!=i){
	numSwap++;
	int temp=detPos(i);
	swap(detPos(i),detPos(temp));
      }
    }
    if (numSwap %2==0)
      return 1;
    else 
      return -1;
}

//assumes swap1 and swap2 are of different spin
complex<double> 
RVBpPsiClass::evaluateRatio_energy(SystemClass &system,int swap1, int swap2)
{
  assert(1==2);
  //  cerr<<"Evaluate ratio energy"<<endl;
  //  mat.SaveInverse();
  ///let's define spin up as swap1
  if (system.x(swap1)!=0){
    swap(swap1,swap2);
  }

  //evaluate the new spin up column
  newCols.resize(1);
  newRows.resize(1);
  newCols[0].resize(u.size());
  newRows[0].resize(u.size());
  colIndices.resize(1);
  colIndices[0]=mat.DetPos[swap1];
  rowIndices.resize(1);
  rowIndices[0]=mat.DetPos[swap2];
  //loop over the spin down particles
  for (int j=0;j<system.x.size();j++){
    if (system.x(j)==1){
      //      u(mat.DetPos(j))=Phi(swap1,j,system);
      newCols[0](mat.DetPos[j])=Phi(swap1,j,system);
    }

  //swap2 has been set to be the spin down value
  //loops over the spin up particles
  //  for (int i=0;i<system.x.size();i++)
    else if (system.x(j)==0){
      //      up(mat.DetPos(i))=Phi(i,swap2,system);
      newRows[0](mat.DetPos[j])=Phi(j,swap2,system);
    }
  }
  complex<double> test_ratio=mat.Ratio_ncol_nrowp(colIndices,rowIndices,newCols,newRows);
  //    complex<double> check_ratio = evaluateRatio_check(system,swap1,swap2);
  //    complex<double> diff=test_ratio-check_ratio;
  //    cerr<<diff<<" "<<check_ratio<<" "<<test_ratio<<endl;
    //  assert((diff*conj(diff)).real()<1e-10);
  //  cerr<<"test ratio is "<<test_ratio<<" "<<endl; //evaluateRatio_check(system,swap1,swap2)<<endl;

  rebuild=false;
  return -test_ratio;
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
    if ( (system.x(i)==-1) || (system.x(i)==2) ){ //HACK!
      upDet++;
      downDet=-1;
      for (int j=0;j<system.x.size();j++){
	if ((system.x(j)==1) || (system.x(j)==2) ){ //HACK!
	  downDet++;
	  //	  cerr<<"GRR GRR"<<upDet<<" "<<downDet<<endl;
	  myMat.M(upDet,downDet)=Phi(i,j,system);
	  myMat.UpPos[i]=upDet;
	  myMat.DownPos[j]=downDet;
	  myMat.DetPos[i]=upDet;
	  myMat.DetPos[j]=downDet;
	}
      }
    }
  }

//   cerr<<"In evaluate ratio "<<endl;
//   for (int i=0;i<mat.UpPos.size();i++){
//     cerr<<i<<" "<<mat.UpPos[i]<<endl;
//   }
//   cerr<<endl;
//   for (int i=0;i<mat.DownPos.size();i++){
//     cerr<<i<<" "<<mat.DownPos[i]<<endl;
//   }

//   cerr<<"Done FILLING NOW"<<endl;

  //  cerr<<"CHECKING NOW"<<endl;
  //  FillDet_check(system,myMat);
  //  cerr<<"done CHECKING NOW"<<endl;
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

