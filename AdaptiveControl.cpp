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
#include <cmath>
#include <fstream>




void rungeKutta(Eigen::VectorXf X0,double t_start, double t_stop) 
{ 
    
    double t_step=0.01;
    long int n = (int)((t_stop - t_start) / t_step); 
    Eigen::VectorXf k1, k2, k3, k4, k5; 
    Eigen::VectorXf k1_dot, k2_dot, k3_dot, k4_dot, k5_dot; 

    Eigen::VectorXf X=X0;
    Eigen::VectorXf dX=Eigen::VectorXf::Zero(6);
    dX(0)=X0(2);
    dX(1)=X0(3);
    X(0)=M_PI/2;
    X(1)=M_PI;


    X(4)=-1;
    X(5)=1;
    Integrator *rk=new Integrator();
    double t=t_start;
    cout<<n;
    
    
    for (int i=1; i<=n; i++) 
    { 
        // Apply Runge Kutta Formulas to find 
        // next value of y

        // dX=rk->diff(X,dX,t); 
        // X=X+ t_step*dX; 
        k1_dot=rk->diff(X,dX,t); 
        k1 = t_step*k1_dot; 
       
        k2_dot=rk->diff(X+0.5*k1,dX+0.5*k1_dot,t+0.5*t_step);
        k2=t_step*k2_dot;
       
            
        k3_dot=rk->diff(X+0.5*k2,dX+0.5*k2_dot,t+0.5*t_step);
        k3=t_step*k3_dot;

        k4_dot=rk->diff(X+k3,dX+k3_dot,t+t_step);
        k4=t_step*k4_dot;
  
        X = X + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        dX=(1.0/6.0)*(k1_dot + 2*k2_dot + 2*k3_dot + k4_dot);
  
        // // Update next value of x 
        if(i%1000==0)
        {
            cout<<"Iteration "<<i<<endl;
            cout<<"delX "<<dX<<endl;
        }
        t = t + t_step; 
        
       
    } 
    vector <pair<double,double>> q1_data=rk->access_q1();
    vector <pair<double,double>> q2_data=rk->access_q2();
    vector <pair<double,double>> q1dot_data=rk->access_q1dot();
    vector <pair<double,double>> q2dot_data=rk->access_q2dot();

    vector <pair<double,double>> m1_data=rk->access_m1();
    vector <pair<double,double>> m2_data=rk->access_m2();

    vector <pair<double,double>> desired_q1_data=rk->access_desired_q1();
    vector <pair<double,double>> desired_q2_data=rk->access_desired_q2();
    vector <pair<double,double>> desired_q1dot_data=rk->access_desired_q1dot();
    vector <pair<double,double>> desired_q2dot_data=rk->access_desired_q2dot();

    ofstream file;
    file.open("q1.csv");

    for(int i=0;i<q1_data.size();i++)
    {
      file<<q1_data[i].second<<","<<q1_data[i].first<<","<<q2_data[i].first<<","<<q1dot_data[i].first
      <<","<<q2dot_data[i].first<<","<<m1_data[i].first<<","<<m2_data[i].first<<","<<desired_q1_data[i].first<<","<<desired_q2_data[i].first<<","
      <<desired_q1dot_data[i].first<<","<<desired_q2dot_data[i].first<<endl;
    }

     
} 
int main()
{
    Eigen::VectorXf X=Eigen::VectorXf::Zero(6);
    Eigen::VectorXf dX=Eigen::VectorXf::Zero(6);
    Eigen::VectorXf new_dX;
    double t_start=0;
    double t_stop=100;

    rungeKutta(X,t_start,t_stop);
    
  
    
}