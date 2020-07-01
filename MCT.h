// a1 - length of link1
// a2 - length of link2
// g  - gravity
// q1 - angle1
// q2 - angle2 
// m1 - mass of link1 at the end
// m2 - mass of link2 at the end

#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <cmath>

class MCT
{
    public:
        MCT(float m1,float m2,float a1,float a2)
        {
            set_m1(m1);
            set_m2(m2);
            set_a1(a1);
            set_a2(a2);
        }

        ~MCT(){}

        void set_m1(float m1) {this->m1=m1;}
        void set_m2(float m2) {this->m2=m2;}
        void set_a1(float a1) {this->a1=a1;}
        void set_a2(float a2) {this->a2=a2;}
        void set_m1_predict(float m1_predict){this->m1=m1_predict;}
        void set_m2_predict(float m2_predict) {this->m2=m2_predict;}
        
        Eigen::MatrixXf M(float q1,float q2)
        {
            Eigen::MatrixXf M(2,2);
            M(0,0)= a1*a1*m1 + (a1*a1 + a2*a2 +2*a1*a2*cos(q2))*m2;
            M(0,1)= a2*(a2 + a1*cos(q2))*m2;
            M(1,0)= m2*a2*(a2 + a1*cos(q2));
            M(1,1)= a2*a2*m2;
            return M;
        }
        Eigen::MatrixXf C(float q1,float q2, float q1dot, float q2dot)
        {
            Eigen::MatrixXf C(2,2);
            C(0,0)=-a1*a2*m2*sin(q2)*q2dot;
            C(0,1)= -a1*a2*m2*sin(q2)*(q1dot +q2dot);
            C(1,0)= a1*a2*m2*sin(q2)*q1dot;
            C(1,1)= 0;
            return C;
        }
        Eigen::MatrixXf T(float q1, float q2)
        {
            Eigen::MatrixXf T(2,1);
            T(0,0)=a1*m1*g*cos(q1) +a1*m2*g*cos(q1) +a2*m2*g*cos(q1+q2);
            T(1,0)=a2*m2*g*cos(q1+q2);
            return T;
        }

    
    private:
        float m1= 3;
        float m2= 2;
        float a1= 2;
        float a2= 3;
        float g= 9.8;
};