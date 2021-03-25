function dadt=Tyson_ode(t,B,a,tau,PP)

%% This is the Tyson 2 variable model ODE

dadt=zeros(size(B));
M=B(1);Pt=B(2);
vm=1;km=0.1;vp=0.5;kp1=10;kp2=0.03;kp3=0.1;Pcrit=0.1;Jp=0.05;

keq=(-100*a*square(2*pi*(1/tau)*t,PP)+100*a)+(200-200*a);

q=2/(1+sqrt(1+8*keq*Pt));
P1=q*Pt;
P2=keq*P1^2;


dadt(1)=vm/(1+(Pt*(1-q)/(2*Pcrit))^2)-km*M;
dadt(2)=vp*M-(kp1*Pt*q+kp2*Pt)/(Jp+Pt)-kp3*Pt;