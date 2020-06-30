
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "MCT.h"
#include "lyap.h"
#include <cmath>

class Integrator
{
public:
    Eigen::VectorXf diff(Eigen::VectorXf X, Eigen::VectorXf dX, float t)
    {
        float q1 = X(0);
        float q2 = X(1);
        float q1dot = X(2);
        float q2dot = X(3);
        float delm1 = X(4);
        float delm2 = X(5);

        float q1ddot = dX(2);
        float q2ddot = dX(3);

        float m1predict = m1 + delm1;
        float m2predict = m2 + delm2;

        cout<<"m1predict "<<m1predict<<endl;
        cout<<"m2predict "<<m2predict<<endl;


        MCT *mctp = new MCT(m1predict, m2predict, a1, a2);
        Eigen::MatrixXf Mp = mctp->M(q1, q2);
        Eigen::MatrixXf Cp = mctp->C(q1, q2, q1dot, q2dot);
        Eigen::MatrixXf Np = mctp->T(q1, q2);
        delete mctp;

        MCT *mct = new MCT(m1, m2, a1, a2);
        Eigen::MatrixXf M = mct->M(q1, q2);
        Eigen::MatrixXf C = mct->C(q1, q2, q1dot, q2dot);
        Eigen::MatrixXf N = mct->T(q1, q2);
        delete mct;


        //calculate tau


        //desired trajectory
        float Dq1 = sin(t);
        float Dq2 = -cos(t);
        float Dqdot1 = cos(t);
        float Dqdot2 = sin(t);
        float Dqddot1 = -sin(t);
        float Dqddot2 = cos(t);
        Eigen::VectorXf DX(6);
        DX(0)=Dq1;
        DX(1)=Dq2;
        DX(2)=Dqdot1;
        DX(3)=Dqdot2;
        DX(4)=Dqddot1;
        DX(5)=Dqddot2;

        Eigen::MatrixXf err(4,1);
        err(0,0)=X(0)-DX(0);
        err(1,0)=X(1)-DX(1);
        err(2,0)=X(2)-DX(2);
        err(3,0)=X(3)-DX(3);
       
    
        
    
        Eigen::MatrixXf W(2, 2);
        W(0, 0) = q1ddot * a1 * a1 + g * cos(q1) * a1;
        W(0, 1) = q1ddot * (a1 * a1 + 2 * cos(q2) * a1 * a2 + a2 * a2) - q2dot * (2 * a1 * a2 * q1dot * sin(q2) + a1 * a2 * q2dot * sin(q2)) + q2ddot * (a2 * a2 + a1 * cos(q2) * a2) + a2 * g * cos(q1 + q2) + a1 * g * cos(q1);
        W(1, 0) = 0;
        W(1, 1) = a2 * (a1 * sin(q2) * q1dot * q1dot + a2 * q1ddot + a2 * q2ddot + g * cos(q1 + q2) + a1 * q1ddot * cos(q2));

        
        //Setting A and Q matrix and lyapunov function
        Eigen::MatrixXf A = Eigen::MatrixXf::Zero(4, 4);
        Eigen::MatrixXf B = Eigen::MatrixXf::Zero(4, 2);
        Eigen::MatrixXf Q = Eigen::MatrixXf::Identity(4, 4);
        Eigen::MatrixXf d_errX(4,1);
        Eigen::MatrixXf del_theta(2,1);

         
        
        A(0, 2) = 1;
        A(1, 3) = 1;
        A(2, 0) = -kp1;
        A(2, 2) = -kd1;
        A(3, 1) = -kp2;
        A(3, 3) = -kd2;


        B(2,0)=1;
        B(3,1)=1;

     

        del_theta(0,0)=delm1;
        del_theta(1,0)=delm2;

        Eigen::MatrixXf Mp_inv=Mp.inverse();
        d_errX=A*err + B*Mp_inv*W*del_theta;
        cout<<"error change "<<d_errX<<endl <<"err"<<err<<endl;

        
        lyap *test = new lyap();
        Eigen::MatrixXf P = test->lyapunov(A.transpose(), Q); 
        

        
        Eigen::MatrixXf mass_ddot=-2*W.transpose()*Mp_inv.transpose()*B.transpose()*P*err;
       
       
        Eigen::VectorXf new_dX(6);
        new_dX(0)=dX(1)+d_errX(0);
        new_dX(1)=dX(2)+d_errX(1);
        new_dX(2)=dX(2)+d_errX(2);
        new_dX(3)=dX(3)+d_errX(3);
        new_dX(4)=mass_ddot(0);
        new_dX(5)=mass_ddot(1);
        
        return new_dX;

    }

private:
    float m1 = 1, m2 = 1, a1 = 1, a2 = 1;
    float g = 9.8;
    float kp1=2,kp2=2, kd1=5,kd2=5;
};
