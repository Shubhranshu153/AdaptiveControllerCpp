#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <cmath>
#include <iostream>

using namespace std;

class lyap
{
    public:
        Eigen::MatrixXf kroneckar(Eigen::MatrixXf I, Eigen::MatrixXf A )
        {
             
            Eigen::MatrixXf P2(I.rows() * A.rows(), I.cols() * A.cols());
            P2.setZero();

            for (int i = 0; i < I.rows(); i++)
            {
                for(int j=0; j<I.cols(); j++)
                    P2.block(i*A.rows(), j*A.cols(), A.rows(), A.cols()) = I(i, j) * A;
            }
           
            return P2;   
        }
        void set_mat_size(int n)
        {
            mat_size=n;
        }
        Eigen::MatrixXf lyapunov(Eigen::MatrixXf A,Eigen::MatrixXf Q_mat)
        {

            set_mat_size(A.rows());
           Eigen::MatrixXf I= Eigen::MatrixXf::Identity (mat_size,mat_size) ;
         
         
     
        
            Eigen::MatrixXf L=kroneckar(I,A);
            Eigen::MatrixXf R=kroneckar(A,I);
            Eigen::MatrixXf T=L+R;


           Eigen::VectorXf Q=Eigen::VectorXf::Zero(mat_size*mat_size);
           
           for(int i=0;i<mat_size;i++)
           {
               for(int j=0;j<mat_size;j++)
               {
                   Q(i*mat_size+j)=Q_mat(j,i);
               }

           }

           Eigen::VectorXf P=-T.inverse()*Q;
          
           
           Eigen::MatrixXf P_mat(mat_size,mat_size);
           for(int i=0;i<mat_size;i++)
           {
               for(int j=0;j<mat_size;j++)
               {
                   P_mat(j,i)=P(i*mat_size+j);
               }

           }
           
            return P_mat;
        }
    private:
        int mat_size; 
};

// int main()
// {
//       Eigen::MatrixXf A= Eigen::MatrixXf::Zero (4,4);
//       Eigen::MatrixXf Q= Eigen::MatrixXf::Identity (4,4);
//            A(0,2)=1;
//            A(1,3)=1;
//            A(2,0)=-2;
//            A(2,2)=-5;
//            A(3,1)=-2;
//            A(3,3)=-5;
//     lyap* test=new lyap();
//      Eigen::MatrixXf P=test->lyapunov(A,Q);

// }