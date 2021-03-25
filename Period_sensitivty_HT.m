function Period_sensitivty_HT
clc;clear all;
%% This is the program to calculate the period sensitivity of Tyson model using HT

dt=0.01;
fs=1/dt;
t=0:dt:2000;
B=0.1*ones(1,2);

dP=0.01;% change in the parameter (10% of the parameter of interest)

[t,A]=ode45(@ode1,t,B,[],dP);

x=A(:,1);
x1=detrend(x);
z = hilbert(x1);
phi=angle(z);
instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
period1=1/mean(instfreq(100000:190000));

period_sensitivity=(period1-24.2)/dP

function dadt=ode1(t,B,dP)

dadt=zeros(size(B));

M=B(1);Pt=B(2);
km=0.1+dP;% Here parameter of interest is km. for other parameter sensitivity, please add...
% dP to that parameter value.

vm=1;vp=0.5;kp1=10;kp2=0.03;kp3=0.1;Pcrit=0.1;Jp=0.05;


keq=200;
q=2/(1+sqrt(1+8*keq*Pt));
P1=q*Pt;
P2=keq*P1^2;


dadt(1)=vm/(1+(Pt*(1-q)/(2*Pcrit))^2)-km*M;
dadt(2)=vp*M-(kp1*Pt*q+kp2*Pt)/(Jp+Pt)-kp3*Pt;