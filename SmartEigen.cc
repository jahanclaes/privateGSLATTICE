#include "SmartEigen.h"

#include "MatrixOps.h"
#include "Blitz.h"
#include <complex>

#include <Eigen/Dense>
#include <vector>
using namespace Eigen;
using namespace std;
complex<double>  SmartEigen::Det()
  {
    return M.determinant();
    //    return MatrixOps::Determinant(M);
  }

void SmartEigen::SaveInverse()
  {
    assert(1==2);
    MInverse_saved=MInverse;
    M_saved=M;
  }


double SmartEigen::GetParity()
{
  vector<int> swapMe;
  for (int i=0;i<UpPos.size();i++){
    if (UpPos[i]!=-1)
      swapMe.push_back(UpPos[i]);
  }

  int countSwaps=0;
  for (int i =swapMe.size()-1;i>=0;i--){
    for (int j=0;j<swapMe.size()-1;j++){
      if (swapMe[j]>swapMe[j+1]){
	swap(swapMe[j],swapMe[j+1]);
	countSwaps++;
      }
    }
  }

  swapMe.clear();
  for (int i=0;i<DownPos.size();i++){
    if (DownPos[i]!=-1)
      swapMe.push_back(DownPos[i]);
  }

  for (int i =swapMe.size()-1;i>=0;i--){
    for (int j=0;j<swapMe.size()-1;j++){
      if (swapMe[j]>swapMe[j+1]){
	swap(swapMe[j],swapMe[j+1]);
	countSwaps++;
      }
    }
  }

  return ((countSwaps % 2) == 0 ? 1 : -1);


}

void SmartEigen::RestoreInverse()
  {
    assert(1==2);
    MInverse=MInverse_saved;
    M=M_saved;
  }

void SmartEigen::Init(int size,int size2)
  {
    M.resize(size,size);
    DetPos.resize(size2);
    UpPos.resize(size2);
    DownPos.resize(size2);
/*     M.resize(size,size); */
/*     MInverse.resize(size,size); */
/*     MInverseu.resize(size); */
/*     MInverse_check.resize(size,size); */
/*     MInverse_saved.resize(size,size); */
/*     MInverse_k.resize(size); */
/*     savedRow.resize(size); */
/*     savedCol.resize(size); */
/*     u.resize(size); */

/*     DetPos.resize(2*size); */
/*     //    cerr<<"DETPOS SIZE IS "<<DetPos.size()<<endl; */
/*     M_saved.resize(size,size); */

/*     //variables for 1 row and 1 col update */
/*     U_r1c1.resize(2,size); */
/*     V_r1c1.resize(size,2); */
/*     Det_UV.resize(2,2); */
/*     MInverseU.resize(2,size); */
/*     MInverseV.resize(size,2); */
/*     VU.resize(size,size); */
/*     MInverseV_DetInverse.resize(size,2); */
  }



  void SmartEigen::CalcAndSaveInverse()
  {
    MInverse=M.inverse();
  }

  complex<double>  SmartEigen::ColRatio(int colIndex,VectorXcd &col)
  {
    complex<double> s = col.dot(MInverse.row(colIndex));
    return s;
  }

  complex<double>  SmartEigen::RowRatio(int colIndex,VectorXcd &col)
  {
    complex<double> s = col.dot(MInverse.col(colIndex));
    return s;
  }

//col is the difference!
  complex<double>  SmartEigen::UpdateColAndInverse(int colIndex,VectorXcd &col)
  {
    VectorXcd u=col-M.col(colIndex);
    M.col(colIndex)=col;
    complex<double> s = col.dot(MInverse.row(colIndex));

    VectorXcd MInverseU = MInverse*u;

    VectorXcd MInverse_k=MInverse.row(colIndex);
    MInverse=MInverse-1.0/s*(MInverseU*MInverse_k.transpose());
    return s;
    /*     u(Range::all())=col(Range::all())-M(colIndex,Range::all());  */
/*     M(colIndex,Range::all())=col(Range::all()); */
/*     firstIndex i; secondIndex j; thirdIndex k; */
/*     MInverse_k(Range::all())=MInverse(Range::all(),colIndex); */
/* /\*     for (int i=0;i<u.size();i++) *\/ */
/* /\*       cerr<<i<<" "<<u(i)<<" "<<MInverse_k(i)<<endl; *\/ */
/*     complex<double>  s=sum(u(i)*MInverse_k(i))+1.0; */
/*     MInverseu=sum(MInverse(j,i)*u(j),j); */
/*     MInverse=MInverse-1.0/s*(MInverseu(j)*MInverse_k(i)); */
/*     return s; */
  }

  complex<double>  SmartEigen::UpdateRowAndInverse(int colIndex,VectorXcd &col)
  {
    VectorXcd u=col-M.row(colIndex).transpose();
    M.row(colIndex)=col;
    complex<double> s = col.dot(MInverse.col(colIndex));
    VectorXcd MInverseU = u.transpose()*MInverse;
    VectorXcd MInverse_k=MInverse.col(colIndex);
    //    cerr<<MInverseU*MInverse_k.transpose()<<endl;
    //MInverse=MInverse-1.0/s*(MInverseU*MInverse_k.transpose());
    MInverse=MInverse-1.0/s*(MInverse_k*MInverseU.transpose());
    return s;

/*     u(Range::all())=col(Range::all())-M(Range::all(),colIndex);  */
/*     M(Range::all(),colIndex)=col(Range::all()); */
/*     firstIndex i; secondIndex j; thirdIndex k; */
/*     MInverse_k(Range::all())=MInverse(colIndex,Range::all()); */
/*     complex<double>  s=sum(u(i)*MInverse_k(i))+1.0; */
/* /\*     for (int i=0;i<u.size();i++) *\/ */
/* /\*       cerr<<i<<" "<<u(i)<<" "<<MInverse_k(i)<<endl; *\/ */
/*     MInverseu=sum(MInverse(i,j)*u(j),j); */
/*     MInverse=MInverse-1.0/s*(MInverseu(i)*MInverse_k(j)); */
/*     return s; */
  }



  ///NEED TO GET THE SIGN CORRECT!
  void SmartEigen::InverseUpdate(vector<int> &colIndices, vector<int> &rowIndices, 
		     vector<blitz::Array<complex<double> ,1> > &newCols,
		     vector<blitz::Array<complex<double> ,1> > &newRows)
  {
/*     int size=M.extent(0); */
/*     int n=colIndices.size(); */
/*     if (n==0) */
/*       assert(1==2); */

/*     U_r1c1.resize(2*n,size); */
/*     V_r1c1.resize(size,2*n); */
/*     Det_UV.resize(2*n,2*n); */
/*     MInverseU.resize(2*n,size); */
/*     MInverseV.resize(size,2*n); */
/*     MInverseV_DetInverse.resize(size,2*n); */
/*     for (int i=0;i<n;i++){ */
/*       V_r1c1(Range::all(),2*i+1)=newRows[i](Range::all())-M(Range::all(),rowIndices[i]); */
/*       V_r1c1(Range::all(),2*i)=0.0; */
/*       V_r1c1(colIndices[i],2*i)=1.0; */
    
/*       U_r1c1(2*i,Range::all())=newCols[i](Range::all())-M(colIndices[i],Range::all()); */

/*       U_r1c1(2*i+1,Range::all())=0.0; */
/*       U_r1c1(2*i+1,rowIndices[i])=1.0; */

/*     } */
/*     for (int i=0;i<colIndices.size();i++)  */
/*       for (int j=0;j<rowIndices.size();j++)  */
/*  	U_r1c1(2*i,rowIndices[j])=0.0;  */



/*     //    MatrixOps::MatrixMultiply(MInverse,V_r1c1,MInverseV); */
/*     MatrixOps::product(MInverse,V_r1c1,MInverseV); */
/*     //    MatrixOps::MatrixMultiply(U_r1c1,MInverse,MInverseU);  */
/*     MatrixOps::product(U_r1c1,MInverse,MInverseU);  */
/*     //    MatrixOps::MatrixMultiply(MInverseU,V_r1c1,Det_UV);  */
/*     MatrixOps::product(MInverseU,V_r1c1,Det_UV);  */
/*     for (int i=0;i<Det_UV.extent(0);i++) */
/*       Det_UV(i,i)=Det_UV(i,i)+1.0;  */
/*     Det_UV=MatrixOps::Inverse(Det_UV); */

/*     //    MatrixOps::MatrixMultiply(MInverseV,Det_UV,MInverseV_DetInverse); */
/*     MatrixOps::product(MInverseV,Det_UV,MInverseV_DetInverse); */
/*     //    MatrixOps::MatrixMultiply(MInverseV_DetInverse,MInverseU,VU); */
/*     MatrixOps::product(MInverseV_DetInverse,MInverseU,VU); */
/*     MInverse=MInverse-VU; */
/*     for (int i=0;i<colIndices.size();i++){ */
/*       int colIndex=colIndices[i]; */
/*       M(colIndex,Range::all())=newCols[i](Range::all()); */
/*     } */

/*     for (int i=0;i<rowIndices.size();i++){ */
/*       int colIndex=rowIndices[i]; */
/*       M(Range::all(),colIndex)=newRows[i](Range::all()); */
/*     } */
/*     //    cerr<<"CHECKING INVERSE:"<<endl; */
/*     ///    CheckInverse(); */
  }

  void SmartEigen::InverseUpdate(vector<int> &colIndices, vector<int> &rowIndices, 
				   Eigen::MatrixXcd &newCols,
				   Eigen::MatrixXcd &newRows)
  {
    int size=M.rows();
    int n=colIndices.size();
     if (n==0) 
       assert(1==2);

     U_r1c1=Eigen::MatrixXcd::Zero(2*n,size);
     V_r1c1=Eigen::MatrixXcd::Zero(size,2*n);

     for (int i=0;i<n;i++){ 
       V_r1c1.col(2*i+1)=newRows.row(i).transpose()-M.col(rowIndices[i]);
       V_r1c1.col(2*i)(colIndices[i])=1.0;
       U_r1c1.row(2*i)=newCols.col(i).transpose()-M.row(colIndices[i]);
       U_r1c1.row(2*i+1)(rowIndices[i])=1.0;
     }
     Eigen::MatrixXcd MInverseV=MInverse*V_r1c1;
     MInverseU=U_r1c1*MInverse;
     Det_UV=MInverseU*V_r1c1;
     Det_UV=Det_UV+Eigen::MatrixXcd::Identity(Det_UV.rows(),Det_UV.cols());
     Eigen::MatrixXcd Det_UVp=Det_UV.inverse();

     //     MInverse=MInverse-MInverse*V_r1c1*Det_UVp*MInverseU;
     MatrixXcd VU=(MInverseV*Det_UVp)*MInverseU;
     MInverse=MInverse-VU;
     for (int i=0;i<colIndices.size();i++)
       M.row(colIndices[i])=newCols.col(i).transpose();

     for (int i=0;i<rowIndices.size();i++)
       M.col(rowIndices[i])=newRows.row(i).transpose();
  }


complex<double> SmartEigen::Ratio_ncol_nrowp(vector<int> &colIndices, 
					     vector<int> &rowIndices,
					     Eigen::MatrixXcd &newCols,
					     Eigen::MatrixXcd &newRows)
{
  int size=M.rows();
  int n=colIndices.size();
  if (n==0)
    return 1.0;

  U_r1c1=Eigen::MatrixXcd::Zero(2*n,size);
  V_r1c1=Eigen::MatrixXcd::Zero(size,2*n);
  for (int i=0;i<n;i++){
    V_r1c1.col(2*i+1)=newRows.row(i).transpose()-M.col(rowIndices[i]);
    V_r1c1.col(2*i)(colIndices[i])=1.0;
    U_r1c1.row(2*i)=newCols.col(i).transpose()-M.row(colIndices[i]);
    U_r1c1.row(2*i+1)(rowIndices[i])=1.0;
  }
  Det_UV.noalias()=U_r1c1*MInverse*V_r1c1;
  Det_UV=Det_UV+Eigen::MatrixXcd::Identity(Det_UV.rows(),Det_UV.cols());
  return Det_UV.determinant();
}





  complex<double>  SmartEigen::Ratio_ncol_nrowp(vector<int> &colIndices, vector<int> &rowIndices, 
			 vector<blitz::Array<complex<double> ,1> > &newCols,
			 vector<blitz::Array<complex<double> ,1> > &newRows)
  {
    assert(1==2);


/*     //    cerr<<"DETPOS SIZE IS "<<DetPos.size()<<endl; */
/*     //    blitz::Array<complex<double> ,2> &U(U_r1c1); */
/*     //    blitz::Array<complex<double> ,2> &V(V_r1c1); */
/*     int size=M.extent(0); */
/*     int n=colIndices.size(); */
/*     if (n==0) */
/*       return 1.0; */
/*     //    cerr<<"SIZES ARE: "<<size<<" "<<2*n<<endl; */
/*     U_r1c1.resize(2*n,size); */
/*     V_r1c1.resize(size,2*n); */
/*     Det_UV.resize(2*n,2*n); */
/*     MInverseU.resize(2*n,size); */
/*     for (int i=0;i<U_r1c1.extent(0);i++) */
/*       for (int j=0;j<U_r1c1.extent(1);j++) */
/*       	U_r1c1(i,j)=0.0; */
/*     for (int i=0;i<V_r1c1.extent(0);i++) */
/*       for (int j=0;j<V_r1c1.extent(1);j++) */
/*        V_r1c1(i,j)=0.0; */
/*     for (int i=0;i<MInverseU.extent(0);i++) */
/*       for (int j=0;j<MInverseU.extent(1);j++) */
/* 	MInverseU(i,j)=0.0; */
/*     for (int i=0;i<Det_UV.extent(0);i++) */
/*       for (int j=0;j<Det_UV.extent(1);j++) */
/*        	Det_UV(i,j)=0.0; */
    

/*     for (int i=0;i<n;i++){ */
/*       V_r1c1(Range::all(),2*i+1)=newRows[i](Range::all())-M(Range::all(),rowIndices[i]); */
/*       V_r1c1(Range::all(),2*i)=0.0; */
/*       V_r1c1(colIndices[i],2*i)=1.0; */
    
/*       for (int j=0;j<U_r1c1.extent(1);j++){ */
/* 	U_r1c1(2*i,j)=newCols[i](j)-M(colIndices[i],j); */
/* 	U_r1c1(2*i+1,j)=0.0; */
/*       }	 */
/*       //      U_r1c1(2*i,Range::all())=newCols[i](Range::all())-M(colIndices[i],Range::all()); */


/*       U_r1c1(2*i+1,rowIndices[i])=1.0; */

/*     } */
/*      for (int i=0;i<colIndices.size();i++)  */
/*        for (int j=0;j<rowIndices.size();j++)  */
/*  	U_r1c1(2*i,rowIndices[j])=0.0;  */

/*     MatrixOps::product(U_r1c1,MInverse,MInverseU);  */
/*     MatrixOps::product(MInverseU,V_r1c1,Det_UV);  */
/*     for (int i=0;i<Det_UV.extent(0);i++) */
/*       Det_UV(i,i)=Det_UV(i,i)+1.0;  */
/*     //    cerr<<Det_UV<<endl;  */
/*     //    cerr<<"DETPOS SIZE IS "<<DetPos.size()<<endl; */
/*     return MatrixOps::Determinant(Det_UV); */
  }


  void SmartEigen::CheckInverse()
  {
    MatrixXcd MInverse_check=M.inverse();
    complex<double> diff=(MInverse_check-MInverse).norm();
    cerr<<"Difference: "<<diff<<endl;
    if (abs(diff)>1e-5){
      cerr<<endl<<endl<<MInverse<<endl<<endl<<endl;
      cerr<<endl<<endl<<MInverse_check<<endl<<endl<<endl;

      assert(1==2);
    }
/*     MInverse_check=MatrixOps::Inverse(M); */
/*     for (int i=0;i<MInverse.extent(0);i++){ */
/*       for (int j=0;j<MInverse.extent(1);j++){ */
/* 	if (!(fabs(abs(MInverse(i,j))-abs(MInverse_check(i,j)))<1e-8)) */
/* 	  cerr<<"INVERSE FAIL: "<<MInverse(i,j)<<" "<<MInverse_check(i,j)<<endl; */
/*       } */
/*       //      cerr<<endl; */
/*     } */
  }
