#ifndef SMARTEIGEN_H
#define SMARTEIGEN_H

#include "MatrixOps.h"
#include "Blitz.h"


#include <Eigen/Dense>
#include <complex>
#include <vector>
//using Eigen::MatrixXcd;
//using Eigen::VectorXcd;


class SmartEigen
{
 public:

  Eigen::MatrixXcd M;
  Eigen::MatrixXcd MInverse;
  Eigen::MatrixXcd MInverse_saved;
  Eigen::MatrixXcd M_saved;



  //  blitz::Array<complex<double> ,2> M;
  //  blitz::Array<complex<double> ,2> MInverse;
  //  blitz::Array<complex<double> ,2> M_saved;

/*   blitz::Array<complex<double> ,2> MInverse_check; */
/*   blitz::Array<complex<double> ,2> MInverse_saved; */
/*   blitz::Array<complex<double> ,1> MInverseu; */
/*   blitz::Array<complex<double> ,1> MInverse_k; */
/*   blitz::Array<complex<double> ,1> savedRow; */
/*   blitz::Array<complex<double> ,1> savedCol; */
/*   blitz::Array<complex<double> ,1> u; */
//  VectorXcd DetPos;
  vector<int> DetPos;
  vector<int> UpPos;
  vector<int> DownPos;
  //  blitz::Array<int,1> DetPos; 

/*   //variables for 1 row and 1 col update */
/*   blitz::Array<complex<double> ,2> U_r1c1; */
/*   blitz::Array<complex<double> ,2> V_r1c1; */
/*   blitz::Array<complex<double> ,2> Det_UV; */
/*   blitz::Array<complex<double> ,2> MInverseU; */
/*   blitz::Array<complex<double> ,2> MInverseV; */
/*   blitz::Array<complex<double> ,2> MInverseV_DetInverse; */
/*   blitz::Array<complex<double> ,2> VU; */
  
  complex<double>  Det();
  void SaveInverse();
  void RestoreInverse();
  void Init(int size);
  void CalcAndSaveInverse();
  complex<double>  ColRatio(int colIndex,Eigen::VectorXcd &col);
  complex<double>  RowRatio(int colIndex,Eigen::VectorXcd &col);
  complex<double>  UpdateColAndInverse(int colIndex,Eigen::VectorXcd &col);
  complex<double>  UpdateRowAndInverse(int colIndex,Eigen::VectorXcd &col);
  void InverseUpdate(vector<int> &colIndices, vector<int> &rowIndices, 
		     vector<blitz::Array<complex<double> ,1> > &newCols,
		     vector<blitz::Array<complex<double> ,1> > &newRows);
  complex<double>  Ratio_ncol_nrowp(vector<int> &colIndices, vector<int> &rowIndices, 
				    vector<blitz::Array<complex<double> ,1> > &newCols,
				    vector<blitz::Array<complex<double> ,1> > &newRows);
  void CheckInverse();

};
#endif
