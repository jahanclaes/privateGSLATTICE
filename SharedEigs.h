#ifndef SHAREDEIGS_H
#define SHAREDEIGS_H
#include "SharedWaveFunctionData.h"
#include "SystemClass.h"
#include "Blitz.h"
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <utility>
#include <fstream>

class SharedEigsClass : public SharedWaveFunctionDataClass
{
 public:
  int NumOrbitals;
  
  void Init(SystemClass &system)
  {
    
    ProcessSystem(system);
  }
  void ProcessSystem(SystemClass &system)
  {
    cerr<<"IN process system"<<" "<<system.rList.size()<<endl;
    GetEigs(system); //HACK FOR HONEYCOMB!
  }
  void SetZero()
  {
    //    eigs=Eigen::MatrixXcd::Zero(eigs.rows(),eigs.cols());
    eigs.setZero();
  }

  Eigen::MatrixXcd eigs;
  void GetEigs(SystemClass &system)
  {
    NumOrbitals=system.kList.size();
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
  

  
};


#endif
