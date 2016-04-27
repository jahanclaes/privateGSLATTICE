#ifndef PROJ_GUTZ_H
#define PROJ_GUTZ_H
#include <iostream>
#include "Blitz.h"
#include "Random/Random.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include "SystemClass.h"
#include "WaveFunction.h"

class ProjectedGutzweilerClass  : public WaveFunctionClass
{
public:
  TinyVector<Array<double,2>,2 > Dets;
  TinyVector<Array<double,2>,2 > MInverse;

  TinyVector<Array<double,1>,2>u;
  Array<double,1> MInverseu;
  Array<double,1> MInverse_k;
  TinyVector<int,6> honeycomb_backup;
  Array<double,2> eigs;
  Array<int,1> DetPos;
  Array<int,1> ParticleOrder;
  void Swap(int i,int j);
  int RealSign;
  int NumSpinUp;
  //currently only the size matters
  void Init(SystemClass &system);
  void UpdateDets(SystemClass &system, int swap1, int swap2);
  void FillDet(SystemClass &system,int spin);
  void GetEigs(SystemClass &system);
  void CheckMInverse();

  complex<double> evaluate(SystemClass &system);
  complex<double> evaluateRatio(SystemClass &system,int swap1, int swap2);
  complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2);
  complex<double> logevaluate(SystemClass &system,int &sign);
  complex<double> evaluateRatio_honeycomb(SystemClass &system,TinyVector<int,6>  &honeycomb_list,
					  TinyVector<int,6> &honeycomb_backup,int amt);
  complex<double> evaluateRatio(SystemClass &system,TinyVector<int,4>  &rhombus_locs,
				TinyVector<int,4> &rhombus_backup);



  void RotateHoneycomb(TinyVector<int,6> &honeycomb,
		       int amt);


};

#endif
