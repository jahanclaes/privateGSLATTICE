#include "SmartPfaffian.h"

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;



double SmartPfaffian::Pfaffian_compute(MatrixXd &m)
{
  MatrixXd mat;
  mat=m;
  int ndim=mat.rows();
  double t1=1.0;
  int kp;
  for (int j=0;j<ndim/2;j++){
    int ndim_r=mat.cols();
    kp=1; double pvt=abs(mat(0,1));
    for (int k=2;k<ndim_r;k++){
      if (abs(mat(0,k))>pvt){
	kp=k;
	pvt=abs(mat(0,k));
      }
    }
    if (kp!=1){
      mat.col(1).swap(mat.col(kp));
      mat.row(1).swap(mat.row(kp));
      t1=-t1;
    }
    t1 = t1 *mat(0,1);
    if (j < ndim/2-1){
      for (int i=2;i<ndim_r;i++){
	if (mat(0,1)!=0.0){
	//	if (abs(mat(0,1))>1e-10){
	  {
	    VectorXd tv=mat.row(1)*mat(i,0)/mat(1,0);
	    mat.row(i)-=tv.transpose();
	  }
	  {	
	    VectorXd tv=mat.col(1)*mat(0,i)/mat(0,1);
	    mat.col(i)-=tv;;
	  }
	}
	else{
	  t1=0.0;
	  break;
	}
      }
    }
    mat=mat.bottomRightCorner(ndim_r-2,ndim_r-2).eval();
  }
  return t1;
}


double SmartPfaffian::Pfaffian()
{
  return  Pfaffian_compute(M);
}


double SmartPfaffian::RowColRatio(int index,Eigen::VectorXd &col)
{
  double total=0.0;
  for (int i=0;i<col.size();i++){
    total+=col(i)*MInverse(index,i);
  }
  return total;
}

void SmartPfaffian::CalcAndStoreInverse()
{
  Eigen::FullPivLU<MatrixXd> lu(M);
  MInverse=lu.inverse(); // M.inverse();
}

// void SmartPfaffian::InverseUpdate(int index, Eigen::VectorXd &col)
// {
//   int size=M.cols();
//   VectorXd V=col.transpose()-M.row(index);
//   cerr<<"A"<<endl;
//   MatrixXd MInverseV=MInverse*V;
//   cerr<<"B"<<endl;
//   MatrixXd Det_UV=V.transpose()*MInverseV;
//   cerr<<"C"<<endl;
//   Det_UV.diagonal().array()+=1.0;
//   cerr<<"D"<<endl;
//   Det_UV=Det_UV.inverse().eval();
//   cerr<<"E"<<endl;
//   MatrixXd VU = MInverseV*Det_UV*V.transpose()*MInverse;
//   cerr<<"F"<<endl;
//   MInverse=MInverse-VU;
//   cerr<<"G"<<endl;
//   M.row(index)=-1*col;
//   cerr<<"H"<<endl;
//   M.col(index)=col;
//   cerr<<"I"<<endl;
//   cerr<<"Actual M Inverse: "<<M.inverse()<<endl;
//   cerr<<"Calculated M Inverse "<<MInverse<<endl;

// }

void SmartPfaffian::InverseUpdate(int index, Eigen::VectorXd &col)
{
  int size=M.cols();
  MatrixXd V(size,2);
  V=MatrixXd::Zero(size,2);
  V.col(0)=-1*col-M.col(index);
  V(index,1)=1.0;
  MatrixXd U(2,size);
  U=MatrixXd::Zero(2,size);
  U.row(1)=1*col.transpose()-M.row(index);
  U(0,index)=1.0;

  MatrixXd MInverseV=MInverse*V;
  MatrixXd Det_UV=U*MInverseV;
  Det_UV.diagonal().array()+=1.0;
  Det_UV=Det_UV.inverse().eval();
  MatrixXd VU = MInverseV*Det_UV*U*MInverse;
  MInverse=MInverse-VU;
  M.row(index)=1*col;
  M.col(index)=-1*col;
}

void SmartPfaffian::CheckInverse()
{
  MatrixXd MInverse_check=M.inverse();
  MatrixXd check=MInverse_check-MInverse;
  cerr<<"The difference is "<<check<<endl;

}
