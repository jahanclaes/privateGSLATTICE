#include "SlaterDet.h"
#include "MatrixOps.h"
// #include "FT.h"
#include <Eigen/Dense>
#include "SmartEigen.h"
#include <algorithm>

/////spin 0 -> first index
////spin 1 -> second index

void 
SlaterDetPsiClass::Init(SystemClass &system,int t_mySpin)
{
  mySpin=t_mySpin;
  NeedFrequentReset=true;
  //  NeedFrequentReset=false;
  SharedEigs.Init(system);
  cerr<<"Starting Init"<<endl;
  NumSpinUp=system.x.size()/2;
  mat.Init(NumSpinUp);
  //  u.resize(NumSpinUp);  
  //  up.resize(NumSpinUp);
  cerr<<"RVB BINS"<<endl;
  //  PairingFunction.Init(system);
  cerr<<"RVB BINS done"<<endl;
  ReadParams=false;
  NumParams=SharedEigs.eigs.cols()*SharedEigs.eigs.rows();
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



//    if (ReadParams){
//      cerr<<"READING IN THE PARMS HERE"<<endl;
//      //    exit(1);
//      ifstream infile;
//      infile.open("params.dat");
//      assert(infile);
//      int num=0;
//      while (!infile.eof()){
//        //       int num;
//        double numGarbage;
//        double val;
//        string garbage;
//        //       infile>>garbage;
//        //       infile>>num;
//        infile>>val;
//        //       infile>>numGarbage;
//        infile>>numGarbage;
//        //HACK!      cerr<<"Setting num val "<<num<<" "<<val<<endl;
//        if (!infile.eof())
//  	PairingFunction.f0[num]=val;
//        num++;
//  	//oldREAD	if (!infile.eof() && num>=2352)
//  	//oldREAD	  PairingFunction.f0[num-2352]=val;
//      }
//     infile.close();

//      for (int i=0;i<system.rList.size();i++)
//        for (int j=0;j<system.rList.size();j++){
//  	int bin=PairingFunction.FindBin(i,j);
//  	cerr<<"READ BIN VALUES: "<<i<<" "<<j<<" "<<PairingFunction.f0[bin].real()<<endl;
//        }


//   }
  
  ReadPairingFunction=false;
//   if (ReadPairingFunction){
//     int bin=0;
//     ifstream infile;
//     infile.open("PairingFunction.dat");
//     assert(infile);
//     while (!infile.eof()){
//       int i;	int j; double tr; double ti;
     
// 	assert(!infile.eof());
// 	infile>>i;
// 	if (!infile.eof()){
// 	assert(!infile.eof());
// 	infile>>j;
// 	bin=PairingFunction.FindBin(i,j);
// 	assert(bin<PairingFunction.f0.size());
// 	assert(!infile.eof());
// 	infile>>PairingFunction.f0[bin].real();
// 	PairingFunction.f0[bin].real(); //=PairingFunction.f0[bin].real()/100.0;
// 	assert(!infile.eof());
// 	infile>>PairingFunction.f0[bin].imag();
// 	PairingFunction.f0[bin].imag(); //=PairingFunction.f0[bin].imag()/100.0;
//       }
//       cerr<<i<<" "<<j<<" "<<PairingFunction.f0[bin]<<" "<<bin<<endl;
//     }
//     infile.close();
//     cerr<<"READ IN PAIRING FUNCTIONS"<<endl;
//   }

//   for (int i=0;i<system.x.size();i++)
//     for (int j=0;j<system.x.size();j++){
//       int bin=PairingFunction.FindBin(i,j);
//       cerr<<"GRR: "<<i<<" "<<j<<" "<<PairingFunction.f0[bin].real()<<" "<<PairingFunction.f0[bin].imag()<<endl;
// 	//<<bin<<endl;
//     }
  cerr<<"done writing "<<endl;
  //  exit(1);
//  NumParams=PairingFunction.f0.size();

//SET NUMPARAMS!!
  evaluate(system);
  cerr<<"Init Done"<<endl;
  //  Array<complex<double>, 1> myDerivs(SharedEigs.NumOrbitals*system.x.size());
  //  AllDerivs(system,myDerivs);

  //  for (int i=0;i<myDerivs.size();i++){
  //    cerr<<i<<" "<<TestDerivs(i,system)<<" "<<myDerivs(i)<<endl;    
  //  }
  //  exit(1);


}


//currently only the size matters
void 
SlaterDetPsiClass::Swap(int i, int j)
{
  swap(mat.DetPos[i],mat.DetPos[j]);
}


void SlaterDetPsiClass::Move(int site, int end_site, int spin)
{
  if (spin==mySpin){
    mat.DetPos[end_site]=mat.DetPos[site];
    mat.DetPos[site]=-1;
  }
    //  if (spin==1){
    //    mat.UpPos[end_site]=mat.UpPos[site];
    //    mat.UpPos[site]=-1;
    //  }
    //  else{
    //    mat.DownPos[end_site]=mat.DownPos[site];
    //    mat.DownPos[site]=-1;
    //  }
}


complex<double>
SlaterDetPsiClass::Deriv(SystemClass &system,int bin)
{
  assert(1==2);
//   complex<double> totalRatio=0.0;
//   for (int i=0;i<system.x.size();i++){
//     for (int j=0;j<system.x.size();j++){
//       if (system.x(i)==1 && system.x(j)==-1){
// 	if (PairingFunction.PhiMap(i,j)==bin){
// 	  totalRatio+=mat.MInverse(mat.DownPos[j],mat.UpPos[i]);
// 	}
// 	else {
	  
// 	}
//       }
//     }
//   }
//  return totalRatio;
}


double SlaterDetPsiClass::GetParam_real(int i)
{
  int numOrbitals=SharedEigs.eigs.rows();
  int numSites=SharedEigs.eigs.cols();
  int orb=i/numSites;
  int site=i % numSites;
  return (SharedEigs.eigs(orb,site)).real();
}

double SlaterDetPsiClass::GetParam_imag(int i)
{
  int numOrbitals=SharedEigs.eigs.rows();
  int numSites=SharedEigs.eigs.cols();
  int orb=i/numSites;
  int site=i % numSites;
  return (SharedEigs.eigs(orb,site)).imag();

}


void SlaterDetPsiClass::SetParam_real(int i, double param)
{
  int numOrbitals=SharedEigs.eigs.rows();
  int numSites=SharedEigs.eigs.cols();
  int orb=i/numSites;
  int site=i % numSites;
  SharedEigs.eigs(orb,site).real()=param;

}
void SlaterDetPsiClass::SetParam_imag(int i, double param)
{
  int numOrbitals=SharedEigs.eigs.rows();
  int numSites=SharedEigs.eigs.cols();
  int orb=i/numSites;
  int site=i % numSites;
  SharedEigs.eigs(orb,site).imag()=param;

}



void
SlaterDetPsiClass::AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs)
{
  AllDerivs(system, derivs,0,derivs.size());

}


void
SlaterDetPsiClass::AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start,int stop)
{
  //  cerr<<"Calling all derivs "<<derivs.size()<<endl;
   int derivPos=start;
   for (int j=0;j<SharedEigs.NumOrbitals;j++){
     for (int site=0;site<system.x.size();site++){
       if (system.x(site)==mySpin || system.x(site)==2){
	 derivs(derivPos)=mat.MInverse(j,mat.DetPos[site]).real();
       }
       else 
	 derivs(derivPos)=0.0;
       derivPos++;
     }
   }
   //   cerr<<"Done with all derivs"<<endl;
}

//   for (int site=0;site<system.x.size();site++){

//     if (system.x(site)==mySpin || system.x(site)==2){
//       Eigen::VectorXcd v=Eigen::VectorXcd::Zero(mat.M.rows());
//       v[mat.DetPos[site]]=1.0;
//       cerr<<"A "<<site<<endl;
//       //      cerr<<mat.MInverse<<endl;
//       double curr=(evaluate_noChange(system)).real();
//       cerr<<"B"<<endl;

//       SharedEigs.eigs(0,site)+=0.01;
//       double post=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(0,site)-=0.02;
//       double pre=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(0,site)+=0.01;
//       cerr<<(post-pre)/0.02/curr<<" "<<mat.MInverse(0,mat.DetPos[site])<<endl;

//       SharedEigs.eigs(1,site)+=0.01;
//       post=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(1,site)-=0.02;
//       pre=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(1,site)+=0.01;
//       cerr<<(post-pre)/0.02/curr<<mat.MInverse(1,mat.DetPos[site])<<endl;


//       SharedEigs.eigs(2,site)+=0.01;
//       post=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(2,site)-=0.02;
//       pre=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(2,site)+=0.01;
//       cerr<<(post-pre)/0.02/curr<<mat.MInverse(2,mat.DetPos[site])<<endl;


//       SharedEigs.eigs(3,site)+=0.01;
//       post=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(3,site)-=0.02;
//       pre=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(3,site)+=0.01;
//       cerr<<(post-pre)/0.02/curr<<mat.MInverse(3,mat.DetPos[site])<<endl;


//       SharedEigs.eigs(4,site)+=0.01;
//       post=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(4,site)-=0.02;
//       pre=(evaluate_noChange(system)).real();
//       SharedEigs.eigs(4,site)+=0.01;
//       cerr<<(post-pre)/0.02/curr<<mat.MInverse(4,mat.DetPos[site])<<endl;




//     }
    


//   }
  
//   assert(1==2);
//   //  mat.CheckInverse();
  
// //   for (int i=start;i<stop;i++)
// //     derivs(i)=0.0;
// //   for (int i=0;i<system.x.size();i++){
// //     for (int j=0;j<system.x.size();j++){
// //             if ((system.x(i)==1 || system.x(i)==2)  && (system.x(j)==-1 || system.x(j)==2)){
// // 	int bin=start+PairingFunction.PhiMap(i,j);
// // 	//	cerr<<"I think the deriv for "<<bin<<" "<<i<<" "<<j<<" "<<system.x(i)<<" "<<system.x(j)<<" "<<" is "<<mat.MInverse(mat.DownPos[j],mat.UpPos[i])<<" "<<TestDerivs(PairingFunction.PhiMap(i,j),system)<<endl;
// // 	if (system.x(i)==2 && system.x(j)==2)
// // 	  derivs(bin)+=mat.MInverse(mat.DownPos[j],mat.UpPos[i]); // +mat.MInverse(mat.DownPos[i],mat.UpPos[j]);
// // 	else
// // 	  derivs(bin)+=mat.MInverse(mat.DownPos[j],mat.UpPos[i]);
// //       }
// //       else {
	
// //       }
// //     }
// //   }

// //   for (int i=0;i<system.x.size();i++){
// //     for (int j=0;j<system.x.size();j++){
// //       if ((system.x(i)==1 || system.x(i)==2)  && (system.x(j)==-1 || system.x(j)==2)){
// // 	int bin=start+PairingFunction.PhiMap(i,j);
// // 	if (abs(derivs(bin)-TestDerivs(PairingFunction.PhiMap(i,j),system))>1e-5){
	  
// // 	  cerr<<"I think the deriv for "<<bin<<" "<<i<<" "<<j<<" "<<system.x(i)<<" "<<system.x(j)<<" "<<" is "<<derivs(bin)<<" "<<TestDerivs(PairingFunction.PhiMap(i,j),system)<<" "<<abs(derivs(bin)-TestDerivs(PairingFunction.PhiMap(i,j),system))<<endl;
// // 	}
// //       }
// //     }
// //   }
// }
double SlaterDetPsiClass::TestDerivs(int derivInt,SystemClass &system)
{
  //  cerr<<"Testing "<<derivInt<<endl;
  double delta=1.0;
  double currParam=GetParam_real(derivInt);

//   for (double myStep=-0.01;myStep<0.01;myStep+=0.001){
//     double step_energy=0.0;
//     double countSteps=0.0;
//     SetParam_real(derivInt,currParam+myStep);
//     cerr<<myStep<<" "<<evaluate_noInverse(system)<<endl;
//   }
  
  SetParam_real(derivInt,currParam+delta);
  double up=evaluate_noChange(system).real();

  SetParam_real(derivInt,currParam-delta);
  double down=evaluate_noChange(system).real();

  SetParam_real(derivInt,currParam);
  double curr=evaluate_noChange(system).real();

  //  cerr<<"Done Testing "<<derivInt<<endl;
  return  (up-down)/(2.0*delta*curr);
}


complex<double>
SlaterDetPsiClass::Phi(int i,int j,SystemClass &system)
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
  //  return PairingFunction.Phi(i,j);
  //  double gutz_correct=PairingFunction.gutz(i,j,system);
  //  if (fabs(gutz_correct-PairingFunction.Phi(i,j).real())>1e-10){
  //    cerr<<"ERROR: "<<gutz_correct<<" "<<PairingFunction.Phi(i,j)<<" "<<i<<" "<<j<<endl;
  //  }
  //  complex<double> toReturn(gutz_correct,0.0);
  //  return toReturn;
}

//void SlaterDetPsiClass::SetParams(double delta3,SystemClass &system)

//{
//  SetParams(0,delta3,system);
//}

//void 
//SlaterDetPsiClass::SetParams(int i,double delta3, SystemClass &system)
//{
//  assert(1==2);
//}


complex<double> 
SlaterDetPsiClass::evaluate_noChange(SystemClass &system)
{
  SmartEigen mat_check;
  mat_check.Init(mat.M.rows());
  FillDet(system,mat_check);
  return mat_check.Det();
}

complex<double> 
SlaterDetPsiClass::evaluate_noInverse(SystemClass &system)
{
  FillDet(system,mat);
  return mat.Det();
}


complex<double>
SlaterDetPsiClass::evaluate(SystemClass &system)
{
  complex<double> myAns=evaluate_noInverse(system);
  mat.CalcAndSaveInverse();
  return myAns;
}


void
SlaterDetPsiClass::Reject(SystemClass &system,int swap1,int swap2)
{
  
}


void SlaterDetPsiClass::Reject(SystemClass &system,int site,int end_site,int spin)
{


}

void 
SlaterDetPsiClass::UpdateDets(SystemClass &system,int swap1, int swap2)
{
  mat.InverseUpdate(colIndices,rowIndices,newCols,newRows);  
}

void 
SlaterDetPsiClass::UpdateDets(SystemClass &system,int site, int end_site,int spin)
{
  //  cerr<<"Starting update det"<<endl;
  if (spin==mySpin){
    //    cerr<<"A"<<endl;
    //    mat.DetPos[end_site]=mat.DetPos[site];
    //    cerr<<"B"<<endl;
    //    mat.DetPos[site]=-1;
    //    cerr<<"C "<<mat.DetPos[end_site]<<endl;
    mat.UpdateRowAndInverse(mat.DetPos[end_site],col);
    //    cerr<<"D"<<endl;
  }
  //  if (spin==1){
    //    //    cerr<<"I'm updating "<<mat.UpPos[end_site]<<endl;
  //    mat.UpdateRowAndInverse(mat.UpPos[end_site],col);
    //    cerr<<"I'm done updating"<<endl;
  //  }
  //  else{
  //    //    cerr<<"I'm d updating "<<mat.DownPos[end_site]<<endl;
  //    mat.UpdateColAndInverse(mat.DownPos[end_site],col);
  //    //    cerr<<"I'm done updating"<<endl;
  //  }
  //  cerr<<"ending  update det"<<endl;

}

complex<double>
SlaterDetPsiClass::evaluateRatio(SystemClass &system,int start, int stop, int spin)
{
  
  //  Array<complex<double> ,1>  derivs;

  //  AllDerivs(system,derivs);
  if (spin!=mySpin)
    return 1.0;
  col.resize(mat.M.rows());
  for (int orb=0;orb<NumSpinUp;orb++){
    col(orb)=SharedEigs.eigs(orb,stop);
  }
  system.Move(stop,start,spin);
  Move(stop,start,spin);
  int pre_sign =Sign(system);
  system.Move(start,stop,spin);
  Move(start, stop, spin);
  int post_sign =Sign(system);
  complex<double> ratio=(mat.RowRatio(mat.DetPos[stop],col))* (double)(pre_sign)*(double)(post_sign);
  //  cerr<<"The ratio is "<<ratio<<" "<<evaluateRatio_check(system,start,stop,spin)<<" "<<pre_sign<<" "<<post_sign<<endl;
  return ratio;
  //  mat.DetPos[stop]=mat.DetPos[start];
  //  mat.DetPos[start]=-1;
  
}

complex<double> 
SlaterDetPsiClass::evaluateRatio_check(SystemClass &system, int site, int end_site,
				     int spin)
{
  SmartEigen mat_check;
  mat_check.Init(mat.M.rows());
  system.Move(end_site,site,spin);
  FillDet(system,mat_check);
  complex<double> pre=mat_check.Det();
  system.Move(site,end_site,spin);
  FillDet(system,mat_check);
  complex<double> post=mat_check.Det();
  return post/pre;
}


//assumes swap1 and swap2 are of different spin
complex<double> 
SlaterDetPsiClass::evaluateRatio(SystemClass &system,int swap1, int swap2)
{

  
  assert(1==2);
  int maxSwap=max(swap1,swap2);
  int minSwap=min(swap1,swap2);
  int mySign=1;
  if ((maxSwap-minSwap) % 2==0)
    mySign=-1;
  else
    mySign=1;

  //  mat.SaveInverse();
  ///let's define spin up as swap1
  if (system.x(swap1)!=0){
    swap(swap1,swap2);
  }

  //evaluate the new spin up column
  newCols.resize(1);
  newRows.resize(1);
  //  newCols[0].resize(u.size());
  //  newRows[0].resize(u.size());
  newCols[0].resize(NumSpinUp);
  newRows[0].resize(NumSpinUp);
  

  colIndices.resize(1);
  colIndices[0]=mat.DetPos[swap1];


  rowIndices.resize(1);
  rowIndices[0]=mat.DetPos[swap2];

  //swap1 has been set to be the spin up value
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
  //  complex<double> check_ratio = evaluateRatio_check(system,swap1,swap2);
//    complex<double> diff=test_ratio-check_ratio;
//    cerr<<diff<<" "<<check_ratio<<" "<<test_ratio<<endl;
      //HACK!      assert((diff*conj(diff)).real()<1e-10);
  //  cerr<<"test ratio is "<<test_ratio<<" "<<endl; //evaluateRatio_check(system,swap1,swap2)<<endl;

  rebuild=false;
  //  cerr<<"CURRENT RATIO IS "<<test_ratio<<endl;
  test_ratio.real()=-1*test_ratio.real(); //*mySign;
  test_ratio.imag()=-1*test_ratio.imag(); //*mySign;
  return test_ratio;
  //return check_ratio;
}

double 
SlaterDetPsiClass::Sign(SystemClass &system)
{

  Array<int,1> detPos(mat.M.rows());
  int ii=-1;
  for (int i=0;i<mat.DetPos.size();i++)
    if (system.x(i)==mySpin || system.x(i)==2){
      ii++;

      //      detPos(ii)=mat.DetPos[i];
      detPos(mat.DetPos[i])=ii;
    }
    int numSwap=0;
    //    for (int k=detPos.size()-1;k>=1;k--){
    for (int step=0;step<detPos.size();step++){
      for (int j=detPos.size()-1;j>=1;j--){
	if (detPos(j)<detPos(j-1)){
	  numSwap++;
	  swap(detPos(j),detPos(j-1));
	}
      }
    }
    if ((numSwap %2)==0)
      return -1;
    else 
      return 1;
}

//assumes swap1 and swap2 are of different spin
complex<double> 
SlaterDetPsiClass::evaluateRatio_energy(SystemClass &system,int swap1, int swap2)
{
  cerr<<"HERE IN EVALUATE RATIO ENERGY"<<endl;
   assert(1==2);
//   //  cerr<<"Evaluate ratio energy"<<endl;
//   //  mat.SaveInverse();
//   ///let's define spin up as swap1
//   if (system.x(swap1)!=0){
//     swap(swap1,swap2);
//   }

//   //evaluate the new spin up column
//   newCols.resize(1);
//   newRows.resize(1);
//   newCols[0].resize(u.size());
//   newRows[0].resize(u.size());
//   colIndices.resize(1);
//   colIndices[0]=mat.DetPos[swap1];
//   rowIndices.resize(1);
//   rowIndices[0]=mat.DetPos[swap2];
//   //loop over the spin down particles
//   for (int j=0;j<system.x.size();j++){
//     if (system.x(j)==1){
//       //      u(mat.DetPos(j))=Phi(swap1,j,system);
//       newCols[0](mat.DetPos[j])=Phi(swap1,j,system);
//     }

//   //swap2 has been set to be the spin down value
//   //loops over the spin up particles
//   //  for (int i=0;i<system.x.size();i++)
//     else if (system.x(j)==0){
//       //      up(mat.DetPos(i))=Phi(i,swap2,system);
//       newRows[0](mat.DetPos[j])=Phi(j,swap2,system);
//     }
//   }
//   complex<double> test_ratio=mat.Ratio_ncol_nrowp(colIndices,rowIndices,newCols,newRows);
//   //    complex<double> check_ratio = evaluateRatio_check(system,swap1,swap2);
//   //    complex<double> diff=test_ratio-check_ratio;
//   //    cerr<<diff<<" "<<check_ratio<<" "<<test_ratio<<endl;
//     //  assert((diff*conj(diff)).real()<1e-10);
//   //  cerr<<"test ratio is "<<test_ratio<<" "<<endl; //evaluateRatio_check(system,swap1,swap2)<<endl;

//   rebuild=false;
//   return -test_ratio;
}



void 
SlaterDetPsiClass::FillDet(SystemClass &system,SmartEigen &myMat)
{
  int detLoc=0;
  for (int ri=0;ri<system.x.size();ri++){
    myMat.DetPos[ri]=-1;
    if (system.x(ri)==mySpin || system.x(ri)==2){
      myMat.DetPos[ri]=detLoc;
      for (int el=0;el<NumSpinUp;el++){
	myMat.M(detLoc,el)=SharedEigs.eigs(el,ri);
      }
      detLoc++;
    }

  }
}





void  
SlaterDetPsiClass::FillDet_check(SystemClass &system,SmartEigen &myMat)
{
  assert(1==2);
//   int upDet=-1;
//   int downDet=-1;
//   bool ok=true;
//   for (int i=0;i<system.x.size();i++){
//     if ((system.x(i)==1) || (system.x(i)==2) ){
//       upDet++;
//       downDet=-1;
//       for (int j=0;j<system.x.size();j++){
// 	if ( (system.x(j)==-1) || (system.x(j)==2)){
// 	  downDet++;
// 	  if (!(myMat.M(myMat.UpPos[i],myMat.DownPos[j])==Phi(i,j,system))){
// 	    cerr<<i<<" "<<j<<" "<<myMat.UpPos[i]<<" "<<myMat.DownPos[j]<<" "<<myMat.M(myMat.UpPos[i],myMat.DownPos[j])<<" "<<Phi(i,j,system)<<" "<<endl;
// 	    ok=false;
// 	  }

// 	}

//       }
//     }
//   }
//   for (int i=0;i<system.x.size();i++)
//     cerr<<system.x(i)<<" ";
//   cerr<<endl;
//   assert(ok);
}
