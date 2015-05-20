#ifndef PAIRINGFUNCTIOMANY_H
#define PAIRINGFUNCTIOMANY_H
#include "SharedWaveFunctionData.h"
#include "SystemClass.h"
#include "Blitz.h"
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <unordered_map>
#include <bitset>
#include <set>
#include <utility>
 
class PairingFunctionMany :public SharedWaveFunctionDataClass
{
 public:

  //  std::tr1::unordered_map<unsigned long, vector<int> > NewMap;
  std::unordered_map<unsigned long, vector<int> > NewMap;
  //  std::tr1::unordered_map<vector<bool>, vector<int> > NewMap2;

  void Init(SystemClass &system)
  {
    ProcessSystem(system);
  }
  void SetParams(double delta3)
  {
  }

  int MaxCorrSize;
  int NumCorrelators;



 complex<double> corr2Val(Array<int,1> &x, int corr)
  {
    int mult=1;
    int bin=0;
    for (int i=0;i<myCorrs[corr].size();i++){
      //      bin=bin+x(myCorrs[corr][i])*mult;
      bin=bin+(x(myCorrs[corr][i])+1)*mult;
      mult*=4;
    }
    //    cerr<<"My value is "<<f0[corr][bin]<<endl;
    //    cerr<<"I am looking at "<<corr<<" "<<bin<<endl;
    assert(corr<f0.size());
    assert(bin<f0[corr].size());
    return f0[corr][bin];
  }


 int corr2Bin(Array<int,1> &x, int corr)
  {
    int mult=1;
    int bin=0;
    for (int i=0;i<myCorrs[corr].size();i++){
      bin=bin+(x(myCorrs[corr][i])+1)*mult;
      mult*=4;
    }
    return bin;
  }


  

  vector<set<int> > correlatorsForSite;
  void BuildCorrelatorsForSites()
  {

    correlatorsForSite.resize(MaxSite+1);
    for (int i=0;i<correlatorsForSite.size();i++){
      assert(i<correlatorsForSite.size());
      correlatorsForSite[i].clear();
    }
    for (int i=0;i<NumCorrelators;i++){
      for (int j=0;j<myCorrs[i].size();j++){
	//	cerr<<"Correlator size is "<<myCorrs[i].size()<<endl;;
	assert(i<myCorrs.size());
	assert(j<myCorrs[i].size());
	assert(myCorrs[i][j]<correlatorsForSite.size());

	correlatorsForSite[myCorrs[i][j]].insert(i);
      }
    }
  }
  


  int MaxSite;
  int NumParams;
  vector<vector<int> > myCorrs;
  vector<pair<int,int> > ParamLoc;
  vector<vector<int> > binLoc;
  void ReadCorrelators(string fileName)
  {
    MaxCorrSize=0;
    MaxSite=0;

    ParamLoc.clear();

    ifstream infile;
    infile.open(fileName.c_str());
    assert(infile);
    infile>>NumCorrelators;
    myCorrs.resize(NumCorrelators);
    binLoc.resize(NumCorrelators);
    NumParams=0;
    for (int i=0;i<NumCorrelators;i++){
      int corrSize;
      infile>>corrSize;
      MaxCorrSize=max(corrSize,MaxCorrSize);

      myCorrs[i].resize(corrSize);

      for (int j=0;j<corrSize;j++){
	infile>>myCorrs[i][j];
	//	NumParams++;
	MaxSite=max(myCorrs[i][j],MaxSite);
	//	binLoc[i][j]=ParamLoc.size();


      }
    }

    f0.resize(NumCorrelators);

    for (int i=0;i<NumCorrelators;i++){
      //      cerr<<"This number is "<<(1<<myCorrs[i].size())<<endl;
      f0[i].resize(1<<(myCorrs[i].size()*2));
      binLoc[i].resize(1<<(myCorrs[i].size()*2));
      for (int j=0;j<binLoc[i].size();j++){
	binLoc[i][j]=NumParams;
	ParamLoc.push_back(make_pair(i,j));
	NumParams++;
      }

    }

    BuildCorrelatorsForSites();


  }
  
  //  vector<complex<double> > f0;
  vector<vector<complex<double> > > f0;

  vector<pair<int,int> > binVals;
  Array<int,2> PhiMap;
  virtual complex<double> Phi(int index, int i,int j)
  {
    //    assert(PhiMap(i,j)>=0 && PhiMap(i,j)<f0.size());
    return f0[PhiMap(i,j)][index];

  }

  void AddBinVal(int i, int j)
  {
    bool found=false;
    if (i>j)
      swap(i,j);
    pair<int,int> myPair(i,j);
    for (int bin=0;bin<binVals.size();bin++){
      if (myPair==binVals[bin])
	found=true;
    }
    if (!found){
      binVals.push_back(myPair);
    }
  }

  int FindBin(int i, int j)
  {
    if (i>j)
      swap(i,j);
    pair<int,int> myPair(i,j);
    for (int bin=0;bin<binVals.size();bin++){
      if (myPair==binVals[bin])
	return bin;
    }
    return -1;
  }

  void ProcessSystem(SystemClass &system)
  {
    ReadCorrelators("corrs.txt");
    for (int i=0;i<NumCorrelators;i++){
      //      if (f0[i].size()==16){ // && i<2){
      for (int j=0;j<f0[i].size();j++){
	f0[i][j]=1.0;
      }

      //if (i==0)
      //	f0[i][0]=1.0;
    	
	//	f0[i][1]=0.0;
	//      }
      //      else if (f0[i].size()==2 && i>=2){
      //	f0[i][0]=0.0;
      //	f0[i][1]=1.0;
      //      }

      //      else {
      //	for (int j=0;j<f0[i].size();j++){
      //	  f0[i][j]=1.0;
      //	}

      for (int j=0;j<f0[i].size();j++){
	f0[i][j]=f0[i][j]+j*0.002+i*0.001;
      }
      
      /*     for (int i=0;i<NumCorrelators;i++) { */
/*       f0[i][0]=1.0; */
/*       f0[i][1]=1.0+0.011*i; */
/*       f0[i][2]=-1+0.009*i; */
/*       f0[i][3]=-1.0+0.01*i; */
/*     } */
    }
/*     for (int i=0;i<NumCorrelators;i++) */
/*       for (int j=0;j<f0[i].size();j++){ */
/* 	f0[i][j]=(i+j)*0.021; */
/*       } */
/*     f0[0][1]=3.0; */
/*         f0[0][2]=0.3; */
/*         f0[0][3]=0.2; */
/*   } */
//    cerr<<"READING THE PARAMETERS"<<endl;
/*     for (int i=0;i<NumCorrelators;i++){ */
/*       for (int j=0;j<f0[i].size();j++){ */
/* 	cerr<<"Parameter: "<<i<<" "<<j<<" "<<f0[i][j]<<endl; */
/*       } */
/*     } */

  }
};


#endif
