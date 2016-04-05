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




// ////////////////////////////////
// ////////////Hubbard//////////
// ////////////////////////////////


void Hubbard::Init(SystemClass &system)
{
  string fileName("U.txt");
  Init(system,fileName);
  
}



void Hubbard::Init(SystemClass &system,string fileName)
{
  ifstream infile;
  infile.open(fileName.c_str());
  if (!infile){
    cerr<<"Could not open hubbard hamiltonian "<<fileName<<endl;
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
  U=1.0;
}


void Hubbard::Set_U(double temp_U)
{
  cerr<<"Setting U to be "<<temp_U<<endl;

  U=temp_U;
}
complex<double> Hubbard::GetEnergyRatio(int site, int end_site, int spin,
					SystemClass &system,
					list<WaveFunctionClass*> &wf_list)
{
  system.Move(site,end_site,spin);
  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    (*wf)->Move(site,end_site,spin);
  complex<double> quick_ratio=1.0;
  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    quick_ratio*=(*wf)->evaluateRatio(system,site,end_site,spin);

  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    (*wf)->Reject(system,site,end_site,spin);
  system.Move(end_site,site,spin);
  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    (*wf)->Move(end_site,site,spin);
  int sign =  ( (system.CountElectrons(site,end_site,spin) % 2) == 0 )  ? 1: -1;
  //  cerr<<"My sign is "<<sign<<endl;
  quick_ratio.real(quick_ratio.real()*sign); // *=sign;
  quick_ratio.imag(quick_ratio.imag()*sign); //*=sign;
  return quick_ratio*-1.0;
}

double Hubbard::Energy(SystemClass &system, 
		       list<WaveFunctionClass*> &wf_list)
{

  complex<double> energy=0.0;
  for (int bond=0;bond<bondList.size();bond++){
    int site=bondList[bond].first;
    int new_site=bondList[bond].second;
    vector<int> spins;

    if (system.x(site)==2){
      spins.push_back(1);
      spins.push_back(-1);
    }
    else if (system.x(site)==0){
      
    }
    else {
      spins.push_back(system.x(site));
    }

    for (int spin_index=0;spin_index<spins.size();spin_index++){

      if ( (system.x(new_site)!=2) && (system.x(new_site)!=spins[spin_index]) ){

	energy+=GetEnergyRatio(site,new_site,spins[spin_index],system,
			       wf_list);
      }

    }

  }
  for (int site=0;site<system.x.size();site++){
    energy+=U*( (system.x(site)==2) ? 1 : 0);
  }
  //  cerr<<"My energy is "<<energy<<endl;
  return energy.real();
}


void Hubbard::AllConnected(list<pair<SpinSwap,double> > &vals, 
			    SystemClass &system, 
			      list<WaveFunctionClass*>  &wf_list)
{
  assert(1==2);

}




// ////////////////////////////////
// ////////////Hopping//////////
// ////////////////////////////////


void Hopping::Init(SystemClass &system)
{
  string fileName("U.txt");
  Init(system,fileName);
  
}



void Hopping::Init(SystemClass &system,string fileName)
{
  ifstream infile;
  infile.open(fileName.c_str());
  if (!infile){
    cerr<<"Could not open hubbard hamiltonian "<<fileName<<endl;
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
  U=1.0;
}


void Hopping::Set_U(double temp_U)
{
  cerr<<"Setting U to be "<<temp_U<<endl;

  U=temp_U;
}
complex<double> Hopping::GetEnergyRatio(int site, int end_site, int spin,
					SystemClass &system,
					list<WaveFunctionClass*> &wf_list)
{
  assert(1==2);
  system.Move(site,end_site,spin);
  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    (*wf)->Move(site,end_site,spin);
  complex<double> quick_ratio=1.0;
  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    quick_ratio*=(*wf)->evaluateRatio(system,site,end_site,spin);

  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    (*wf)->Reject(system,site,end_site,spin);
  system.Move(end_site,site,spin);
  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    (*wf)->Move(end_site,site,spin);

  int sign =  ( (system.CountElectrons(site,end_site,spin) % 2) == 0 )  ? 1: -1;
  quick_ratio.real()*=sign;
  quick_ratio.imag()*=sign;
  return quick_ratio;
}

double Hopping::Energy(SystemClass &system, 
		       list<WaveFunctionClass*> &wf_list)
{

  complex<double> energy=0.0;
  for (int bond=0;bond<bondList.size();bond++){
    int site=bondList[bond].first;
    int new_site=bondList[bond].second;
    vector<int> spins;

    if (system.x(site)==2){
      spins.push_back(1);
      spins.push_back(-1);
    }
    else if (system.x(site)==0){
      
    }
    else {
      spins.push_back(system.x(site));
    }

    for (int spin_index=0;spin_index<spins.size();spin_index++){

      if ( (system.x(new_site)!=2) && (system.x(new_site)!=spins[spin_index]) ){

	energy+=GetEnergyRatio(site,new_site,spins[spin_index],system,
			       wf_list);
      }

    }

  }
  for (int site=0;site<system.x.size();site++){
    energy+=U*( (system.x(site)==2) ? 1 : 0);
  }
  //  cerr<<"My energy is "<<energy<<endl;
  return energy.real();
}


void Hopping::AllConnected(list<pair<SpinSwap,double> > &vals, 
			    SystemClass &system, 
			      list<WaveFunctionClass*>  &wf_list)
{
  assert(1==2);

}




// ////////////////////////////////
// ////////////RingExchange//////////
// ////////////////////////////////


void RingExchange::Init(SystemClass &system)
{
  string fileName("U.txt");
  Init(system,fileName);
  
}



void RingExchange::Init(SystemClass &system,string fileName)
{
  ifstream infile;
  infile.open(fileName.c_str());
  if (!infile){
    cerr<<"Could not open hubbard hamiltonian "<<fileName<<endl;
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
  U=1.0;
}


void RingExchange::Set_U(double temp_U)
{
  cerr<<"Setting U to be "<<temp_U<<endl;

  U=temp_U;
}
complex<double> RingExchange::GetEnergyRatio(int site, int end_site, int spin,
					SystemClass &system,
					list<WaveFunctionClass*> &wf_list)
{
  system.Move(site,end_site,spin);
  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    (*wf)->Move(site,end_site,spin);
  complex<double> quick_ratio=1.0;
  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    quick_ratio*=(*wf)->evaluateRatio(system,site,end_site,spin);

  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    (*wf)->Reject(system,site,end_site,spin);
  system.Move(end_site,site,spin);
  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
    (*wf)->Move(end_site,site,spin);
  return quick_ratio;
}

double RingExchange::Energy(SystemClass &system, 
		       list<WaveFunctionClass*> &wf_list)
{

  complex<double> energy=0.0;
  for (int bond=0;bond<bondList.size();bond++){
    int site=bondList[bond].first;
    int new_site=bondList[bond].second;
    vector<int> spins;

    if (system.x(site)==2){
      spins.push_back(1);
      spins.push_back(-1);
    }
    else if (system.x(site)==0){
      
    }
    else {
      spins.push_back(system.x(site));
    }

    for (int spin_index=0;spin_index<spins.size();spin_index++){

      if ( (system.x(new_site)!=2) && (system.x(new_site)!=spins[spin_index]) ){

	energy+=GetEnergyRatio(site,new_site,spins[spin_index],system,
			       wf_list);
      }

    }

  }
  for (int site=0;site<system.x.size();site++){
    energy+=U*( (system.x(site)==2) ? 1 : 0);
  }
  //  cerr<<"My energy is "<<energy<<endl;
  return energy.real();
}


void RingExchange::AllConnected(list<pair<SpinSwap,double> > &vals, 
			    SystemClass &system, 
			      list<WaveFunctionClass*>  &wf_list)
{
  assert(1==2);

}


