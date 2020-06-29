
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "MCT.h"
#include <cmath>

class integrator
{
    public:
        Eigen::MatrixXf diff(Eigen::VectorXf X, Eigen::VectorXf dX)
        {
            float q1=X(0);
            float q2=X(1);
            float q1dot=X(2);
            float q2dot=X(3);
            float delm1=X(4);
            float delm2=X(5);

            float m1predict=m1 +delm1;
            float m2predict=m2 +delm2;

            MCT *mctp =new MCT(m1predict,m2predict,a1,a2);
            Eigen::MatrixXf Mp=mctp->M(q1,q2);
            Eigen::MatrixXf Cp=mctp->C(q1,q2,q1dot,q2dot);
            Eigen::MatrixXf Np=mctp->T(q1,q2);
            delete mctp;

            MCT *mct =new MCT(m1,m2,a1,a2);
            Eigen::MatrixXf M=mct->M(q1,q2);
            Eigen::MatrixXf C=mct->C(q1,q2,q1dot,q2dot);
            Eigen::MatrixXf N=mct->T(q1,q2);
            delete mct;




        }
    private:
        float m1=1,m2=1,a1=1,a2=1;


}