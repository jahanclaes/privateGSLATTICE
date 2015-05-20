#include <iostream>
#include "Communication/Communication.h"
#include "Random/Random.h"
#include "OptimizeBothClass.h"
#include <omp.h>
#include <bitset>
#include "VMC.h"
#include "DMC.h"
#include "input.h"





int main(int argc, char **argv)
{
  cerr<<"Run beginning"<<endl;
  COMM::Init(argc, argv);
  int th_id, nthreads;
#pragma omp parallel private(th_id) shared(nthreads)
  {
    th_id = omp_get_thread_num();
    #pragma omp critical
    {
      cout << "Hello World from thread " << th_id << '\n';
    }
    #pragma omp barrier
 
    #pragma omp master
    {
      nthreads = omp_get_num_threads();
      cout << "There are " << nthreads << " threads" << '\n';
    }
  }

  string fileName(argv[1]);
  cerr<<"File name is "<<fileName<<endl;
  ifstream infile;
  infile.open(fileName.c_str());
  InputClass input;
  input.Read(infile);
  infile.close();

  assert(input.IsVariable("RunType"));
  string myRunType=input.GetVariable("RunType");
  if (myRunType=="OPT"){
    cerr<<"Running opt"<<endl;
    VMCDriverClass VMCDriver;
    VMCDriver.RunMultipleOpt(input);

  }
  else {
    VMCDriverClass VMCDriver;
    VMCDriver.RunSingleVMC(input);
    //    VMCDriver.RunMultipleOpt(input);
    assert(1==2);
  }
  exit(1);
//   infile.open(fileName.c_str());

  
//   string line;
//   infile>>line;
//   cerr<<"line is "<<line<<endl;
//   if (line=="OPT"){
//     //    cerr<<"Running opt"<<endl;
//     //    VMCDriverClass VMCDriver;
//     //    VMCDriver.RunMultipleOpt(infile);
//     //    VMCDriver.RunFiniteT_exact(infile);
//   }
//   else {
//     VMCDriverClass VMCDriver;
//     VMCDriver.RunSingleVMC(infile);

//   }
  
  //  BuildHamiltonian();
  //  DMCDriverClass DMCDriver;
  //  DMCDriver.RunMultipleDMC_CT(J2);
  //  RunSingleVMC(J2);
  //  RunMultipleVMC(J2);
  //  RunSingleVMC(J2);
  exit(1);
}
