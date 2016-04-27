#ifndef SMARTPFAFFIAN_H
#define SMARTPFAFFIAN_H

#include <Eigen/Dense>
#include <iostream>
#include <vector>
using namespace std;

#include <complex>
class SmartPfaffian
{
  
 public:
  Eigen::MatrixXd M;
  Eigen::MatrixXd MInverse;
  Eigen::MatrixXd M_saved;
  Eigen::MatrixXd MInverse_saved;
  vector<int> DetPos;
  
   void Init(int size)
   {
     M.resize(size,size);
     MInverse.resize(size,size);
     M_saved.resize(size,size);
     MInverse_saved.resize(size,size);
     DetPos.resize(2*size);
     for (int i=0;i<DetPos.size();i++)
       DetPos[i]=-1;
   }



/*   Array<complex<double> ,2> M; */
/*   Array<complex<double> ,2> M_saved; */
/*   Array<complex<double> ,2> MInverse; */
/*   Array<complex<double> ,2> MInverse_check; */
/*   Array<complex<double> ,2> MInverse_saved; */
/*   Array<complex<double> ,1> MInverseu; */
/*   Array<complex<double> ,1> MInverse_k; */
/*   Array<complex<double> ,1> savedRow; */
/*   Array<complex<double> ,1> savedCol; */
/*   Array<complex<double> ,1> u; */
/*   Array<int,1> DetPos; */

/*   //variables for 1 row and 1 col update */
/*   Array<complex<double> ,2> U_r1c1; */
/*   Array<complex<double> ,2> V_r1c1; */
/*   Array<complex<double> ,2> Det_UV; */
/*   Array<complex<double> ,2> MInverseU; */
/*   Array<complex<double> ,2> MInverseV; */
/*   Array<complex<double> ,2> MInverseV_DetInverse; */
/*   Array<complex<double> ,2> VU; */
 

  double Pfaffian();
  double Pfaffian_compute(Eigen::MatrixXd &m);

  void SaveInverse() 
  { 
     MInverse_saved=MInverse; 
     M_saved=M; 
   }


   void RestoreInverse() 
   { 
     MInverse=MInverse_saved; 
     M=M_saved; 
   } 

   void CalcAndStoreInverse();
   
   
   double RowColRatio(int rowIndex,Eigen::VectorXd &col);
   

/*   void Init(int size) */
/*   { */
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
/*   } */


   void InverseUpdate(int index, Eigen::VectorXd &col);


/*    ///NEED TO GET THE SIGN CORRECT!  */
/*    void InverseUpdate(vector<int> &colIndices, vector<int> &rowIndices,   */
/*  		     vector<Array<complex<double> ,1> > &newCols,  */
/*  		     vector<Array<complex<double> ,1> > &newRows)  */
/*    {  */
/*      int size=M.extent(0);  */
/*      int n=colIndices.size();  */
/*      if (n==0)  */
/*        assert(1==2);  */

/*      U_r1c1.resize(2*n,size);  */
/*      V_r1c1.resize(size,2*n);  */
/*      Det_UV.resize(2*n,2*n);  */
/*      MInverseU.resize(2*n,size);  */
/*      MInverseV.resize(size,2*n);  */
/*      MInverseV_DetInverse.resize(size,2*n);  */
/*      for (int i=0;i<n;i++){  */
/*        V_r1c1(Range::all(),2*i+1)=newRows[i](Range::all())-M(Range::all(),rowIndices[i]);  */
/*        V_r1c1(Range::all(),2*i)=0.0;  */
/*        V_r1c1(colIndices[i],2*i)=1.0;  */
    
/*        U_r1c1(2*i,Range::all())=newCols[i](Range::all())-M(colIndices[i],Range::all());  */

/*        U_r1c1(2*i+1,Range::all())=0.0;  */
/*        U_r1c1(2*i+1,rowIndices[i])=1.0;  */

/*      }  */
/*      for (int i=0;i<colIndices.size();i++)   */
/*        for (int j=0;j<rowIndices.size();j++)   */
/*   	U_r1c1(2*i,rowIndices[j])=0.0;   */



/*      //    MatrixOps::MatrixMultiply(MInverse,V_r1c1,MInverseV);  */
/*      MatrixOps::product(MInverse,V_r1c1,MInverseV);  */
/*      //    MatrixOps::MatrixMultiply(U_r1c1,MInverse,MInverseU);   */
/*      MatrixOps::product(U_r1c1,MInverse,MInverseU);   */
/*      //    MatrixOps::MatrixMultiply(MInverseU,V_r1c1,Det_UV);   */
/*      MatrixOps::product(MInverseU,V_r1c1,Det_UV);   */
/*      for (int i=0;i<Det_UV.extent(0);i++)  */
/*        Det_UV(i,i)=Det_UV(i,i)+1.0;   */
/*      Det_UV=MatrixOps::Inverse(Det_UV);  */

/*      //    MatrixOps::MatrixMultiply(MInverseV,Det_UV,MInverseV_DetInverse);  */
/*      MatrixOps::product(MInverseV,Det_UV,MInverseV_DetInverse);  */
/*      //    MatrixOps::MatrixMultiply(MInverseV_DetInverse,MInverseU,VU);  */
/*      MatrixOps::product(MInverseV_DetInverse,MInverseU,VU);  */
/*      MInverse=MInverse-VU;  */
/*      for (int i=0;i<colIndices.size();i++){  */
/*        int colIndex=colIndices[i];  */
/*        M(colIndex,Range::all())=newCols[i](Range::all());  */
/*      }  */

/*      for (int i=0;i<rowIndices.size();i++){  */
/*        int colIndex=rowIndices[i];  */
/*        M(Range::all(),colIndex)=newRows[i](Range::all());  */
/*      }  */
/*      //    cerr<<"CHECKING INVERSE:"<<endl;  */
/*      ///    CheckInverse();  */
/*    }  */






/*   ///NEED TO GET THE SIGN CORRECT! */

/*   complex<double>  Ratio_ncol_nrow(vector<int> &colIndices, vector<int> &rowIndices,  */
/* 			 vector<Array<complex<double> ,1> > &newCols, */
/* 			 vector<Array<complex<double> ,1> > &newRows) */
/*   { */
/*     //    cerr<<"DETPOS SIZE IS "<<DetPos.size()<<endl; */
/*     //    Array<complex<double> ,2> &U(U_r1c1); */
/*     //    Array<complex<double> ,2> &V(V_r1c1); */
/*     int size=M.extent(0); */
/*     int n=colIndices.size(); */
/*     if (n==0) */
/*       return 1.0; */
/*     //    cerr<<"SIZES ARE: "<<size<<" "<<2*n<<endl; */
/*     U_r1c1.resize(2*n,size); */
/*     V_r1c1.resize(size,2*n); */
/*     Det_UV.resize(2*n,2*n); */
/*     MInverseU.resize(2*n,size); */
/*     for (int i=0;i<n;i++){ */
/*       V_r1c1(Range::all(),2*i+1)=newRows[i](Range::all())-M(Range::all(),rowIndices[i]); */
/*       V_r1c1(Range::all(),2*i)=0.0; */
/*       V_r1c1(colIndices[i],2*i)=1.0; */
    
/*       U_r1c1(2*i,Range::all())=newCols[i](Range::all())-M(colIndices[i],Range::all()); */

/*       U_r1c1(2*i+1,Range::all())=0.0; */
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
/*   } */


/*   complex<double>  Ratio_ncol_nrowp(vector<int> &colIndices, vector<int> &rowIndices,  */
/* 			 vector<Array<complex<double> ,1> > &newCols, */
/* 			 vector<Array<complex<double> ,1> > &newRows) */
/*   { */
/*     //    cerr<<"DETPOS SIZE IS "<<DetPos.size()<<endl; */
/*     //    Array<complex<double> ,2> &U(U_r1c1); */
/*     //    Array<complex<double> ,2> &V(V_r1c1); */
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
/*   } */


/*   complex<double>  UpdateColAndInverse(int colIndex,Array<complex<double> ,1> &col) */
/*   { */
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
/*   } */

/*   complex<double>  UpdateRowAndInverse(int colIndex,Array<complex<double> ,1> &col) */
/*   { */
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
/*   } */
/*   complex<double>  RowRatio(int colIndex,Array<complex<double> ,1> &col) */
/*   { */
/*     firstIndex i; secondIndex j; thirdIndex k; */
/*     MInverse_k(Range::all())=MInverse(colIndex,Range::all()); */
/*     complex<double>  s=sum(col(i)*MInverse_k(i)); */
/*     return s; */
/*   } */
   void CheckInverse();

/*   void CheckInverse() */
/*   { */
/*     MInverse_check=MatrixOps::Inverse(M); */
/*     for (int i=0;i<MInverse.extent(0);i++){ */
/*       for (int j=0;j<MInverse.extent(1);j++){ */
/* 	if (!(fabs(abs(MInverse(i,j))-abs(MInverse_check(i,j)))<1e-8)) */
/* 	  cerr<<"INVERSE FAIL: "<<MInverse(i,j)<<" "<<MInverse_check(i,j)<<endl; */
/*       } */
/*       //      cerr<<endl; */
/*     } */
/*   } */

};
#endif
