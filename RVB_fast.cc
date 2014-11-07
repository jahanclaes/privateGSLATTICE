#include "RVB_fast.h"
#include "MatrixOps.h"
// #include "FT.h"
#include "SmartMatrix.h"
#include <algorithm>

/////spin 0 -> first index
////spin 1 -> second index

void 
RVBFastPsiClass::Init(SystemClass &system)
{
    NeedFrequentReset=true;
  //  NeedFrequentReset=false;
  cerr<<"Starting Init"<<endl;
  NumSpinUp=system.x.size()/2;
  mat.Init(NumSpinUp);
  u.resize(NumSpinUp);  
  up.resize(NumSpinUp);
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
     cerr<<"READING IN THE PARMS"<<endl;
     //    exit(1);
     ifstream infile;
     infile.open("params.dat");
     assert(infile);
     while (!infile.eof()){
       int num;
       double numGarbage;
       double val;
       string garbage;
       infile>>garbage;
       infile>>num;
       infile>>val;
       infile>>numGarbage;
       infile>>numGarbage;
       //HACK!      cerr<<"Setting num val "<<num<<" "<<val<<endl;
       if (!infile.eof())
 	PairingFunction.f0[num]=val;
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
  
  ReadPairingFunction=true;
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
	//	double junk;
	//	infile>>junk;
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
RVBFastPsiClass::Swap(int i, int j)
{
  swap(mat.DetPos(i),mat.DetPos(j));
}



complex<double>
RVBFastPsiClass::Deriv(SystemClass &system,int bin)
{
  complex<double> totalRatio=0.0;
  for (int i=0;i<system.x.size();i++){
    for (int j=0;j<system.x.size();j++){
      if (system.x(i)==0 && system.x(j)==1){
	if (PairingFunction.PhiMap(i,j)==bin){
	  totalRatio+=mat.MInverse(mat.DetPos(j),mat.DetPos(i));
	}
	else {
	  
	}
      }
    }
  }
  return totalRatio;
}


double RVBFastPsiClass::GetParam_real(int i)
{
  return PairingFunction.f0[i].real();

}

double RVBFastPsiClass::GetParam_imag(int i)
{
  return PairingFunction.f0[i].imag();

}


void RVBFastPsiClass::SetParam_real(int i, double param)
{
  PairingFunction.f0[i].real()=param;

}
void RVBFastPsiClass::SetParam_imag(int i, double param)
{
  PairingFunction.f0[i].imag()=param;

}




void
RVBFastPsiClass::AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs)
{
  AllDerivs(system, derivs,0,derivs.size());

}

void
RVBFastPsiClass::AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start,int stop)
{


  for (int i=start;i<stop;i++)
    derivs(i)=0.0;
  for (int i=0;i<system.x.size();i++){
    for (int j=0;j<system.x.size();j++){
      if (system.x(i)==0 && system.x(j)==1){
	int bin=start+PairingFunction.PhiMap(i,j);
	derivs(bin)+=mat.MInverse(mat.DetPos(j),mat.DetPos(i));

      }
      else {
	
      }
    }
  }

}


complex<double>
RVBFastPsiClass::Phi(int i,int j,SystemClass &system)
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
  return PairingFunction.Phi(i,j);
  double gutz_correct=PairingFunction.gutz(i,j,system);
  if (fabs(gutz_correct-PairingFunction.Phi(i,j).real())>1e-10){
    cerr<<"ERROR: "<<gutz_correct<<" "<<PairingFunction.Phi(i,j)<<" "<<i<<" "<<j<<" "<<gutz_correct/PairingFunction.Phi(i,j).real()<<endl;
  }
  complex<double> toReturn(gutz_correct,0.0);
  return toReturn;
}

void RVBFastPsiClass::SetParams(double delta3,SystemClass &system)

{
  SetParams(0,delta3,system);
}

void 
RVBFastPsiClass::SetParams(int i,double delta3, SystemClass &system)
{
  assert(1==2);
}

complex<double> 
RVBFastPsiClass::evaluate_noInverse(SystemClass &system)
{
  FillDet(system,mat);
  return mat.Det();
}


complex<double>
RVBFastPsiClass::evaluate(SystemClass &system)
{
  complex<double> myAns=evaluate_noInverse(system);
  mat.CalcAndSaveInverse();
  return myAns;
}


void
RVBFastPsiClass::Reject(SystemClass &system,int swap1,int swap2)
{

}

void 
RVBFastPsiClass::UpdateDets(SystemClass &system,int swap1, int swap2)
{
  mat.InverseUpdate(colIndices,rowIndices,newCols,newRows);  
}


//assumes swap1 and swap2 are of different spin
complex<double> 
RVBFastPsiClass::evaluateRatio(SystemClass &system,int swap1, int swap2)
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
  if (system.x(swap1)!=0){
    swap(swap1,swap2);
  }

  //evaluate the new spin up column
  newCols.resize(1);
  newRows.resize(1);
  newCols[0].resize(u.size());
  newRows[0].resize(u.size());

  colIndices.resize(1);
  colIndices[0]=mat.DetPos(swap1);


  rowIndices.resize(1);
  rowIndices[0]=mat.DetPos(swap2);

  //swap1 has been set to be the spin up value
  //loop over the spin down particles
  for (int j=0;j<system.x.size();j++){
    if (system.x(j)==1){
      //      u(mat.DetPos(j))=Phi(swap1,j,system);
      newCols[0](mat.DetPos(j))=Phi(swap1,j,system);
    }

  //swap2 has been set to be the spin down value
  //loops over the spin up particles
  //  for (int i=0;i<system.x.size();i++)
    else if (system.x(j)==0){
      //      up(mat.DetPos(i))=Phi(i,swap2,system);
      newRows[0](mat.DetPos(j))=Phi(j,swap2,system);
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
RVBFastPsiClass::Sign(SystemClass &system)
{
    Array<int,1> detPos(mat.DetPos.size());
    for (int i=0;i<mat.DetPos.size();i++)
      if (system.x(i)==0)
	detPos(i)=mat.DetPos(i);
      else 
	detPos(i)=mat.DetPos(i)+system.x.size()/2;
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
RVBFastPsiClass::evaluateRatio_energy(SystemClass &system,int swap1, int swap2)
{
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
  colIndices[0]=mat.DetPos(swap1);
  rowIndices.resize(1);
  rowIndices[0]=mat.DetPos(swap2);
  //loop over the spin down particles
  for (int j=0;j<system.x.size();j++){
    if (system.x(j)==1){
      //      u(mat.DetPos(j))=Phi(swap1,j,system);
      newCols[0](mat.DetPos(j))=Phi(swap1,j,system);
    }

  //swap2 has been set to be the spin down value
  //loops over the spin up particles
  //  for (int i=0;i<system.x.size();i++)
    else if (system.x(j)==0){
      //      up(mat.DetPos(i))=Phi(i,swap2,system);
      newRows[0](mat.DetPos(j))=Phi(j,swap2,system);
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
RVBFastPsiClass::FillDet(SystemClass &system,SmartMatrix &myMat)
{
  int upDet=-1;
  int downDet=-1;
  for (int i=0;i<system.x.size();i++){
    if (system.x(i)==0){
      upDet++;
      downDet=-1;
      for (int j=0;j<system.x.size();j++){
	if (system.x(j)==1){
	  downDet++;
	  myMat.M(upDet,downDet)=Phi(i,j,system);
	  myMat.DetPos(i)=upDet;
	  myMat.DetPos(j)=downDet;
	}
      }
    }
  }
}


