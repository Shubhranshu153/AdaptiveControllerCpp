% For q1 vs q1des
global A B q1 q2 q1dot1 q2dot1 delm1 delm2 q1dot2 q2dot2 t_old;
global kp1 kp2 kd1 kd2 m1 m2 a1 a2 g m1hat m2hat time1 j  V Vdot ;
j=1;
q1=pi/2;
q2=pi;
q1dot1=0;
q2dot1=0;

m1hat=0;
m2hat=0;

q1dot2=0;
q2dot2=0;

%model parameters CAn be changed as desired
m1=3;
m2=2;
a1=2;
a2=3;
g=9.81;
global m1predict;global m2predict;
m1predict=2;
m2predict=3;
 

% choosing parameters for q1
kp1 =2;
kd1=10;
kp2 =2;
kd2=10;


A=[0 0 1 0; 0 0 0 1; -kp1 0 -kd1 0; 0 -kp2 0 -kd2];
B=[0 0; 0 0; 1 0; 0 1];

delt=0.1;
time_elasped=0;
i=1;
q1array =zeros(100);

time_span=0:0.1:100;

delm1=m1predict-m1;
delm2=m2predict-m2;

 i=1;
while i<2
    
      xinit =[q1 q2 q1dot1 q2dot1 delm1 delm2];  
      [t,x]=ode45(@ODE_45_Func,time_span,xinit);
 
      
      figure(1);
      hold on
      plot(t, x(:,1),'r');
      plot(t, sin(t),'b');
      title('q1 and q1 desired vs time');
      xlabel('t');
      ylabel('q1');
      legend('q1','q1 desired');
      
      figure(2);
      hold on
      plot( x(:,2),'r');
      plot( -cos(t),'b');
      title('q2 and q2 desired vs time');
      xlabel('t');
      ylabel('q2');
      legend('q2','q2 desired');
      
      figure(3);
      hold on
      plot(t, x(:,3),'r');
      plot(t, cos(t),'b');
      title('q1dot and q1dot desired vs time');
      xlabel('t');
      ylabel('q1dot');
      legend('q1dot','q1dot desired');
      
      figure(4);
      hold on
      plot(t, x(:,4),'r');
      plot(t, sin(t),'b');
      title('q2dot and q2dot desired vs time');
      xlabel('t');
      ylabel('q2dot');
      legend('q2dot','q2dot desired');
      
      figure(5);
      hold on
      plot( m1hat,'r');
      plot( m2hat,'b');
      title('m1predict m2predict vs time');
      xlabel('t');
      ylabel('mass');
      legend('m1predict','m2predict');
      
      
      figure(6);
      hold on
      plot(  V,'r');
      title('V vs time');
      xlabel('t');
      ylabel('V');
      
    
      figure(7);
      hold on 
      plot(Vdot,'r');
      title('Vdot vs time');
      xlabel('t');
      ylabel('Vdot');
      

      
     i=i+1; 
end
plot(q1array);
