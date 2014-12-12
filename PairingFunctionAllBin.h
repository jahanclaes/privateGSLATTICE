#ifndef PAIRINGFUNCTIONALLBIN_H
#define PAIRINGFUNCTIONALLBIN_H
#include "SharedWaveFunctionData.h"
#include "SystemClass.h"
#include "Blitz.h"
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <utility>
#include <fstream>

class PairingFunctionAllBin : public SharedWaveFunctionDataClass
{
 public:
  void Init(SystemClass &system)
  {

    ProcessSystem(system);
  }
  vector<complex<double> > f0;
  vector<pair<int,int> > binVals;
  Array<int,2> PhiMap;
  virtual complex<double> Phi(int i,int j)
  {
    //    cerr<<i<<" "<<j<<" "<<PhiMap(i,j)<<" "<<f0.size()<<endl;
    //    assert(PhiMap(i,j)>=0 && PhiMap(i,j)<f0.size());
    return f0[PhiMap(i,j)];

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
      //      cerr<<"Adding val: "<<i<<" "<<j<<endl;
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
    cerr<<"IN process system"<<" "<<system.rList.size()<<endl;
    binVals.clear();
    
    for (int j=0;j<system.rList.size();j++){  
      for (int i=0;i<system.rList.size();i++){  
	AddBinVal(i,j);
      }
    }
    cerr<<"There are "<<binVals.size()<<" "<<"bins"<<endl;
    cerr<<"system size is "<<system.rList.size()<<endl;
    PhiMap.resize(system.rList.size(),system.rList.size());
    for (int i=0;i<system.rList.size();i++)
      for (int j=0;j<system.rList.size();j++){
	int bin=FindBin(i,j);
	assert(bin!=-1);
	assert(bin>=0);
	
	PhiMap(i,j)=bin;
      }
    //    for (int i=0;i<binVals.size();i++)
    //      cerr<<binVals[i].first<<" "<<binVals[i].second<<endl;
    f0.resize(binVals.size());    
    cerr<<"Number of bins are "<<f0.size()<<endl;
    //HACK FOR INITIALIZING!
    for (int i=0;i<f0.size();i++){
      f0[i]=1+i*0.001;
    }

    GetEigs(system); //HACK FOR HONEYCOMB!
    gutz_build(system); //HACK FOR HONEYCOMB!
  }


  Array<complex<double> ,2> eigs;
  void GetEigs(SystemClass &system)
  {
    eigs.resize(system.kList.size(),system.rList.size());
    ifstream infile;
    infile.open("eigs.txt");
    assert(infile);
    cerr<<"Klist size is "<<system.kList.size()<<endl;
    cerr<<"rlist size is "<<system.rList.size()<<endl;

    for (int i=0;i<system.kList.size();i++)
      for (int j=0;j<system.rList.size();j++){
	infile>>eigs(i,j).real();
	//	infile>>eigs(i,j).imag();
	eigs(i,j).imag()=0;
      }
    infile.close();
  }
  
  
  double gutz(int i,int j,SystemClass &system)
  {
    
    int N=system.x.size();
    //    assert(i!=j);
    double total=0.0;
    for (int ki=0;ki<eigs.extent(0);ki++){
      assert(ki<eigs.extent(0));
      assert(i<eigs.extent(1));
      assert(j<eigs.extent(1));
      //      total+=eigs(ki,j)*eigs(ki,i);
      total+=eigs(ki,j).real()*eigs(ki,i).real()+eigs(ki,j).imag()*eigs(ki,i).imag();
    }
    return total;
  }

  void gutz_build(SystemClass &system)
  {


    f0.resize(binVals.size());
    for (int i=0;i<f0.size();i++)
      f0[i]=0.0;
    for (int i=0;i<system.rList.size();i++){
      for (int j=0;j<system.rList.size();j++){
	int bin=FindBin(i,j);
	assert(bin!=-1);
	f0[bin]=gutz(i,j,system);
	//HACK! 
	//	f0[bin].imag()=1.9;
	//	f0[bin].real()+=0.3;
      }
    }
  }

  
};


#endif
