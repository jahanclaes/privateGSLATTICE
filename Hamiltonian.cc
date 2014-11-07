#include "Hamiltonian.h"


#include "SystemClass.h"
#include "WaveFunction.h"
#include <fstream>

void HamiltonianClass::Set_J(double temp_J)
{
   assert(1==2);
}



double HamiltonianClass::Energy(SystemClass &system, 
				list<WaveFunctionClass*> &wf_list)
{

  assert(1==2);
}

void HamiltonianClass::Set_J2(double t_J2)
{
  cerr<<"NOT SETTING J2"<<endl;

}

void HamiltonianClass::Init(SystemClass &system)
{
  assert(1==2);
}



void HamiltonianClass::AllConnected(list<pair<SpinSwap,double> > &vals, 
				    SystemClass &system, 
				    list<WaveFunctionClass*>  &wf_list)
{

  assert(1==2);
}



// ////////////////////////////////
// ////////////Heisenberg//////////
// ////////////////////////////////


void Heisenberg::Init(SystemClass &system)
{
  string fileName("J.txt");
  Init(system,fileName);
  
}



void Heisenberg::Init(SystemClass &system,string fileName)
{
  ifstream infile;
  infile.open(fileName.c_str());
  if (!infile){
    cerr<<"Could not open heisenberg hamiltonian "<<fileName<<endl;
    cerr<<"Aborting!"<<endl;
    exit(1);
  }
  
  while (!infile.eof()){
    int i;
    int j;
    infile>>i;
    infile>>j;
    if (!infile.eof())
      bondList.push_back(make_pair(i,j));
  }
  infile.close();
  J=1.0;
}


void Heisenberg::Set_J(double temp_J)
{
  cerr<<"Setting J to be "<<temp_J<<endl;

  J=temp_J;
}


double Heisenberg::Energy(SystemClass &system, 
			  list<WaveFunctionClass*> &wf_list)
{
  double Jz=J; // 1.0; HACK! HACK! NEED TO SWITCH FOR PYROCHLORE
  double Jx=J; //1.0;
  complex<double> total_Jz=0.0;
  complex<double> total_Jx=0.0;
  complex<double> total=0.0;
  for (int counter=0;counter<bondList.size();counter++){
    int i=bondList[counter].first;
    int j=bondList[counter].second;
    if (system.x(i)==system.x(j)){
      total_Jz+=1.0*Jz;
    }
    else {
      total_Jz-=1.0*Jz;
      system.Swap(i,j);
      for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	(*wf)->Swap(i,j);
      complex<double> quick_ratio=1.0;
      for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	quick_ratio*=(*wf)->evaluateRatio(system,i,j);
      for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	(*wf)->Reject(system,i,j);
      system.Swap(i,j);
      for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	(*wf)->Swap(i,j);
      //HACK!      cerr<<"My quick ratio is "<<quick_ratio<<endl;
      total_Jx=total_Jx+2.0 * quick_ratio * Jx;
    }
  }
  //The total energy eventually needs to be real so let's just some
  //the real part.  For a check it would be nice to know that the
  //complex part was really zero.

  //the 4.0 for the hbar/2 piece
  //  cerr<<"Computing energy for "<<Name<<endl;
  //  cerr<<"with a J of "<<J<<" "<<J*total.real()/4.0<<endl;
  //  cerr<<"Energies are "<<total_Jz.real()<<" "<<total_Jx.real()<<endl;
  term1+=total_Jz.real()/4.;
  term2+=total_Jx.real()/4.;
  //cerr<<"ENERGIES: "<<total_Jz.real()/4.<<" "<<total_Jx.real()/4.0<<endl;
  return (total_Jz+total_Jx).real()/4.0;
  //  return total.real()/4.0;// *J
}


void Heisenberg::AllConnected(list<pair<SpinSwap,double> > &vals, 
			    SystemClass &system, 
			      list<WaveFunctionClass*>  &wf_list)
{
  assert(1==2);

}


