function dxdt=ODE_45_Func(t,x)
global A B q1 q2 q1dot1 q2dot1 delm1 delm2 q1dot2 q2dot2 t_old; 
global kp1 kp2 kd1 kd2 m1 m2 a1 a2 g time1 j m1hat m2hat  V Vdot ;

global m1predict;
global m2predict;


time1(j)=t; 
dxdt = zeros(6,1);



W =[q1dot2*(a1^2)+g*a1*cos(q1),q1dot2*(a1^2+2*cos(q2)*a1*a2 +a2^2) - q2dot1*(2*a1*a2*q1dot1*sin(q2)+a1*a2*q2dot1*sin(q2))+q2dot2*(a2^2+a1*cos(q2)*a2)+a2*g*cos(q1+q2)+a1*g*cos(q1);0,a2*(a1*sin(q2)*q1dot1^2 +a2*q1dot2 +a2*q2dot2 +g*cos(q1+q2)+a1*q1dot2*cos(q2))];  
Mp=[a1^2*m1predict + (a1^2+a2^2+2*a1*a2*cos(q2))*m2predict, a2*(a2+a1*cos(q2))*m2predict; m2predict*a2*(a2 +a1*cos(q2)), a2^2*m2predict];

qd1= sin(t);
qd2 = -cos(t);
Dqd11=cos(t);
Dqd21=sin(t);
Dqd12= -sin(t);
Dqd22= cos(t);

M_hat_inv=inv(Mp);
Q=eye(4);
A1 = [0 0 1 0;0 0 0 1; -kp1 0 -kd1 0; 0 -kp2 0 -kd2];
P1 = lyap(A1',Q);
DX=[qd1;qd2;Dqd11;Dqd21];

del_dxdt = zeros(4,1);
del_dxdt(1:4)=A1*(x(1:4)-DX(1:4)) +B*(M_hat_inv*W)*[delm1;delm2];
dxdt(1)=Dqd11+del_dxdt(1);
dxdt(2)=Dqd21 + del_dxdt(2);
dxdt(3)=Dqd12 + del_dxdt(3);
dxdt(4)=Dqd22 + del_dxdt(4);



k2= -2*W'*M_hat_inv*(B')*P1*(x(1:4)-DX(1:4));

time1(j)=t;
display(t);

dxdt(5) = k2(1,1);
dxdt(6) = k2(2,1);
    
delm1 =x(5);
delm2 =x(6);   
q1 =x(1); 
q2 = x(2);
q1dot1 =x(3);
q2dot1=x(4);
q1dot2=dxdt(3);
q2dot2=dxdt(4);

m1predict=m1+delm1;
m2predict=m2+delm2;

m1hat(j) =m1predict;
m2hat(j)=m2predict;
j=j+1;
delX =x(1:4)-DX(1:4);
deltheta= [delm1;delm2];
V(j)=delX'*P1*delX + 0.5*(deltheta')*deltheta;
Vdot(j) = -delX'*(Q)*delX;

display(x(5));
display(x(6));


end