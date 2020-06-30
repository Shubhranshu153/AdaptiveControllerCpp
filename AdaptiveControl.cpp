// a1 - length of link1
// a2 - length of link2
// g  - gravity
// q1 - angle1
// q2 - angle2 
// m1 - mass of link1 at the end
// m2 - mass of link2 at the end

#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <iostream>
#include "integrator.h"
void rungeKutta(Eigen::VectorXf X0,float t_start, float t_stop) 
{ 
    
    float t_step=0.01;
    int n = (int)((t_stop - t_start) / t_step); 
    Eigen::VectorXf k1, k2, k3, k4, k5; 
    Eigen::VectorXf k1_dot, k2_dot, k3_dot, k4_dot, k5_dot; 

    Eigen::VectorXf X=X0;
    Eigen::VectorXf dX=Eigen::VectorXf::Zero(6);
    dX(0)=X0(2);
    dX(1)=X0(3);
    Integrator *rk=new Integrator();
    float t=t_start;
    cout<<n;
    
    for (int i=1; i<=10; i++) 
    { 
        // Apply Runge Kutta Formulas to find 
        // next value of y
        k1_dot=rk->diff(X,dX,t); 
        k1 = t_step*k1_dot; 
       
        k2_dot=rk->diff(X+0.5*k1,dX+0.5*k1_dot,t+0.5*t_step);
        k2=t_step*k2_dot;
       
            
        k3_dot=rk->diff(X+0.5*k2,dX+0.5*k2_dot,t+0.5*t_step);
        k3=t_step*k3_dot;

        k4_dot=rk->diff(X+k3,dX+k3_dot,t+t_step);
        k4=t_step*k4_dot;
  
        X = X + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        dX=dX + (1.0/6.0)*(k1_dot + 2*k2_dot + 2*k3_dot + k4_dot);
  
        // // Update next value of x 
        t = t + t_step; 
        
        cout<<"change of delm1 "<<X(4)<<"\nchange of delm2 "<<X(5)<<endl;
    } 
  
     
} 
int main()
{
    Eigen::VectorXf X=Eigen::VectorXf::Zero(6);
    Eigen::VectorXf dX=Eigen::VectorXf::Zero(6);
    Eigen::VectorXf new_dX;
    float t_start=0;
    float t_stop=10;

    rungeKutta(X,t_start,t_stop);
    
  
    
}