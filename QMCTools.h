#ifndef QMC_TOOLS_H
#define QMC_TOOLS_H
class QMCTools
{
 public:
static double CalcEnergy(vector<OptimizeBothClass*> &VMC_vector)
{
  double totalE=0.0;
  for (int i=0;i<VMC_vector.size();i++){
    for (list<HamiltonianClass*>::iterator ham_iter=VMC_vector[i]->Ham.begin();
	 ham_iter!=VMC_vector[i]->Ham.end();ham_iter++){
      double vecE=0.0;
      vecE=(*ham_iter)->Energy(VMC_vector[i]->System,VMC_vector[i]->wf_list);
      totalE+=vecE;
    }
  }
  cerr<<"The current energy is "<<totalE/VMC_vector.size()<<endl;
  return totalE/VMC_vector.size();
}


//Now we need to decide where we are putting the reconfigured values.
static void Branch(vector<OptimizeBothClass*> &VMC_vector,
	    vector<OptimizeBothClass*> &VMC_vector_branch,
	    vector<double> &weight,
	    RandomClass &Random)
{

  assert(VMC_vector.size()==VMC_vector_branch.size());
  for (int i=0;i<VMC_vector.size();i++)
    VMC_vector_branch[i]->Copy(*VMC_vector[i]);


//   for (int i=0;i<VMC_vector.size();i++)
//     VMC_vector[i]->Copy(*VMC_vector_branch[i]);
//   return; 

 
  vector<double> cumulant;
  //    double w=0.0;
  double w=0.0;
  //    for (list<pair<SpinSwap,double> >::iterator iter=vals.begin();iter!=vals.end();iter++){
  for (int i=0;i<weight.size();i++){
    assert(weight[i]>-1e-25);
    w+=weight[i];
    cumulant.push_back(w);
  }
  cerr<<"My total weight is "<<w<<endl;
  //really could just multiply waht we are lookign for by w
  for (int i=0;i<cumulant.size();i++)
    cumulant[i]=cumulant[i]/w;



  
  for (int j=0;j<VMC_vector.size();j++){
    double toFind=Random.ranf();
    int i=0;
    while (toFind>cumulant[i])
      i++;
    VMC_vector[j]->Copy(*VMC_vector_branch[i]);
    //    VMC_vector[j]->Copy(*VMC_vector_branch[(j+1) % VMC_vector.size()]);
  }
	
}


static int Sample(list<pair<SpinSwap,double> > &vals,
	   RandomClass &Random,
	   double &w)
{
  vector<double> cumulant;
  //    double w=0.0;
  w=0.0;
  for (list<pair<SpinSwap,double> >::iterator iter=vals.begin();iter!=vals.end();iter++){
    //      cerr<<"My vales are "<<(*iter).second<<" "<<(*iter).first.spin1<<" "<<(*iter).first.spin2<<endl;
    if (((*iter).second<-1e-10))
      cerr<<"My vales are "<<(*iter).second<<" "<<(*iter).first.spin1<<" "<<(*iter).first.spin2<<endl;
    assert((*iter).second>-1e-10);
    w+=(*iter).second;
    cumulant.push_back(w);
  }
  //really could just multiply waht we are lookign for by w
  for (int i=0;i<cumulant.size();i++)
    cumulant[i]=cumulant[i]/w;
  double toFind=Random.ranf();
  int i=0;
  while (toFind>cumulant[i])
    i++;
  return i;
}
};


#endif
