
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "MCT.h"
#include "lyap.h"
#include <cmath>
#include <vector>


class Integrator
{
public:
    Eigen::VectorXf diff(Eigen::VectorXf X, Eigen::VectorXf dX, double t)
    {
        static int i=1;
        double q1 = X(0);

        double q2 = X(1);
        double q1dot = X(2);
        double q2dot = X(3);
        double delm1 = X(4);
        double delm2 = X(5);

        double q1ddot = dX(2);
        double q2ddot = dX(3);

        double m1predict = m1 + delm1;
        double m2predict = m2 + delm2;
        if(i%1000==0)
            cout<<"Mass "<<m1predict<<" "<<m2predict<<endl;\
        i++;

        MCT *mctp = new MCT(m1predict, m2predict, a1, a2);
        Eigen::MatrixXf Mp = mctp->M(q1, q2);
        delete mctp;

        //desired trajectory
        double Dq1 = sin(t);
        double Dq2 = -cos(t);
        double Dqdot1 = cos(t);
        double Dqdot2 = sin(t);
        double Dqddot1 = -sin(t);
        double Dqddot2 = cos(t);
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


        B(2,0)=1;
        B(3,1)=1;

     

        del_theta(0,0)=delm1;
        del_theta(1,0)=delm2;
      //  cout<<"\n\nA "<<A<<endl;

         
        
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
      //  cout<<"\n\nA "<<A<<endl;

        Eigen::MatrixXf Mp_inv=Mp.inverse();
        d_errX=A*err + B*Mp_inv*W*del_theta;
       

        
        lyap *test = new lyap();
        Eigen::MatrixXf P = test->lyapunov(A.transpose(), Q); 
     //   cout<<"P "<<P<<endl;

        
        Eigen::MatrixXf mass_ddot=-2*W.transpose()*Mp_inv*B.transpose()*P*err;
   
        Eigen::VectorXf new_dX(6);
        
       // cout<<"mass change "<<mass_ddot<<endl
        new_dX(0)=DX(2)+d_errX(0);
        new_dX(1)=DX(3)+d_errX(1);
        new_dX(2)=DX(4)+d_errX(2);
        new_dX(3)=DX(5)+d_errX(3);
        new_dX(4)=mass_ddot(0);
        new_dX(5)=mass_ddot(1);

        q1_plt.push_back(std::pair<double,double>(X(0),t));
        q2_plt.push_back(std::pair<double,double>(X(1),t));
        
        
        return new_dX;

    }

    vector<std::pair<double,double>> access_q1()
    {
        return q1_plt;
    }

    vector<std::pair<double,double>> access_q2()
    {
        return q2_plt;
    }


private:
    double m1 = 3, m2 = 2, a1 = 2, a2 = 3;
    double g = 9.8;
    double kp1=2,kp2=2, kd1=10,kd2=10;
    vector <std::pair<double,double>> q1_plt;
    vector <std::pair<double,double>> q2_plt;
};

