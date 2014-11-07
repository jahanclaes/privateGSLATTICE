#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <list>
#include <vector>
#include "Blitz.h"

class SystemClass;
class WaveFunctionClass;

class SpinSwap
{
 public:
  int spin1;
  int spin2;
};


class HamiltonianClass
{
 public:
  double term1;
  double term2;
  string Name;
  Array<double,2> EnergyHist;
  double NumHistCounts;
  virtual double Energy(SystemClass &system,
			list<WaveFunctionClass*> &wf_list);
  virtual void Init(SystemClass &system);
  HamiltonianClass()
    {
    }
  HamiltonianClass(string myName)
    {
      Name=myName;
    }


  virtual void Set_J2(double t_J2);
  virtual void Set_J(double temp_J);

  virtual void AllConnected(list<pair<SpinSwap,double> > &vals,
			    SystemClass &system,
			    list<WaveFunctionClass*>  &wf_list);


};


class Heisenberg : public HamiltonianClass
{
 public:
 Heisenberg(string myName): HamiltonianClass(myName)
  {
  }

  double J;
  double Energy(SystemClass &system,
		list<WaveFunctionClass*> &wf_list);
  vector<pair<int,int> > bondList;
  void Init(SystemClass &system);
  void Init(SystemClass &system,string fileName);
  void Set_J(double temp_J);
  void AllConnected(list<pair<SpinSwap,double> > &vals,
		    SystemClass &system,
		    list<WaveFunctionClass*>  &wf_list);

};


#endif
