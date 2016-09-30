#ifndef PAIRINGFUNCTIOMANY_H
#define PAIRINGFUNCTIOMANY_H
#include "SharedWaveFunctionData.h"
#include "SystemClass.h"
#include "Blitz.h"
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <tr1/unordered_map>
#include <bitset>
#include <set>
#include <cmath>
#include <utility>
 
class PairingFunctionMany :public SharedWaveFunctionDataClass
{
 public:

  std::tr1::unordered_map<unsigned long, vector<int> > NewMap;
  //  std::tr1::unordered_map<vector<bool>, vector<int> > NewMap2;

  void Init(SystemClass &system)
  {
    ProcessSystem(system);
  }
  void SetParams(double delta3)
  {
  }

  int NumCorrsCover;
  int MaxCorrSize;
  int NumCorrelators;



 complex<double> corr2Val(Array<int,1> &x, int corr)
  {
    int mult=1;
    int bin=0;
    for (int i=0;i<myCorrs[corr].size();i++){
      bin=bin+x(myCorrs[corr][i])*mult;
      mult*=2;
    }
    return f0[corr][bin];
  }


 int corr2Bin(Array<int,1> &x, int corr)
  {
    int mult=1;
    int bin=0;
    for (int i=0;i<myCorrs[corr].size();i++){
      bin=bin+x(myCorrs[corr][i])*mult;
      mult*=2;
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
    infile>>NumCorrsCover;
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
      f0[i].resize(1<<myCorrs[i].size());
      binLoc[i].resize(1<<myCorrs[i].size());
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
      for (int j=0;j<f0[i].size();j++){
        if (i>=NumCorrsCover)
	        f0[i][j]=1;
        else
            f0[i][j]=pow(-1,i+j);
      }
    }
    for (int i=0;i<NumCorrelators;i++){
      for (int j=0;j<f0[i].size();j++){
	    cerr<<"Parameter: "<<i<<" "<<j<<" "<<f0[i][j]<<endl;
      }
    }
  }
};


#endif
