function Tyson_Phaseslip
clc;clear all;
%% This is the program to find the phase slip wrt to external signal.

dt=0.01;
Fs=1/dt;
t_end=6000;
t1=0:dt:t_end;
B=0*ones(1,2);

%% Input, amplitude and period of forcing signal
T_E=23;% forcing period
a=0.19;% forcing amplitude
%% **********************************************

[t,A]=ode23s(@ode1,t1,B,[],T_E,a);

x=A(:,1);
%% *******************************
x2=x;
x1=detrend(x2);  
z = hilbert(x1);
phi=(angle(z));
% period calculation
x=phi;
N=0;
zerIdx=[];
for i=1:length(x)-1
   if ((x(i)>N && x(i+1)<N))
      zerIdx(end+1)=i; % save index of zero-crossing
   end
end
t1=t(zerIdx); 
for j=1:length(t1)-1
    period1(j)=t1(j+1)-t1(j);
end
period1;
M1=mean(period1);%find out the mean period

%% period calculation (cycle by cycle)
cycle=round(t_end/T_E);
period_bin=[];
% 
for td=1:cycle-20
    ynm=period1(td)*ones(1,(T_E/0.01));
   period_bin=[period_bin ynm];   
end

period_bin=period_bin';

slip(1)=0;


for i=1:length(period_bin);
    slip(i+1)=slip(i)+dt*((period_bin(i)-T_E)./period_bin(i));
end

  plot(t(1:length(slip))/24,slip/24);
  xlabel('day')
  ylabel('slip(days)');
  slip_per_day=(period1-T_E)*24/period1
 
function dadt=ode1(t,B,T_E,a);

dadt=zeros(size(B));
M=B(1);Pt=B(2);
vm=1;km=0.1;vp=0.5;kp1=10;kp2=0.03;kp3=0.1;Pcrit=0.1;Jp=0.05;

keq=(-100*a*square(2*pi*(1/T_E)*t,50)+100*a)+(200-200*a);
 
q=2/(1+sqrt(1+8*keq*Pt));

dadt(1)=vm/(1+(Pt*(1-q)/(2*Pcrit))^2)-km*M;
dadt(2)=vp*M-(kp1*Pt*q+kp2*Pt)/(Jp+Pt)-kp3*Pt;