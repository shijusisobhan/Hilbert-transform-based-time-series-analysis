function PRC_Tyson_HT
clc;clear all;
%% This is the program to construct the PRC of Tyson et al model using HT
time=0:0.01:1500;

number_of_points_required=12;%How many ponits required for PRC

N_rep=1;%Number of PRC cycle required

Number_cycle=5;%Number of cycle after which phase calculating

T_ref=3.5 ;%specify the circadian time of peak value of the fist varible

P_light=198; %Intensity of perturbation(For type 1 change P_light=185) 

N_light=0; %Unpertubed intensity

P_time=2;%perturbation duration (For type 1 change P_time=0.7)

B=0.1*ones(1,2);
%% 

option = odeset('RelTol', 1e-8);

[t1,A1]=ode23s(@ode2,time,B,option,[],N_light,N_light,[]);

%******************************************************************************
%% Period claculation

x1=find(time>500);
x=A1(:,1);
A_ref=x;
x_new=x(x1);
x=x_new;
N=mean(x_new);
zerIdx=[];
for i=1:length(x)-1
   if ((x(i)<N && x(i+1)>N))
      zerIdx(end+1)=i; % save index of zero-crossing
   end
end

t1=time(zerIdx); % display timepoints of zero-crossings

for j=1:length(t1)-1
    period1(j)=t1(j+1)-t1(j);
end
period1;
period=mean(period1)%find out the mean period

kk1=find(time>=1000 & time<=1000+period);

Time_new=time(kk1);

kk2=find((A_ref(kk1))==max(A_ref(kk1)));

CT0=Time_new(kk2)-T_ref*(period/24);

T_initial=CT0;


%% *********************************************************************************

time_interval=(period*N_rep/number_of_points_required);
Tp=T_initial:time_interval:T_initial+period*N_rep; %Tp is the time at which perturbation applied

for i=1:length(Tp)
  perturbation_time=Tp(i);

[time,A]=ode23s(@ode2,time,B,option,perturbation_time,P_light,N_light,P_time);

T_phase_calculation=T_initial+Number_cycle*period;%Phase calculation after N cycle

%%  phase difference calculation using hilbert transform
T1=find(time>=T_phase_calculation);%phase shift calculate after 20 cycle
X=A1(:,2);Y=A(:,2);
x_pt=X(T1);y_pt=Y(T1);
T_new=time(T1);

h1=hilbert(detrend(x_pt));
h2=hilbert(detrend(y_pt));

p1 = unwrap(angle(h1)); p2 =unwrap(angle(h2)); % Instantaneous phase (wrapped)
tau=period/(2*pi);

p=mean(p2-p1);

 phase_shift(i)=p*tau;
 
 if phase_shift(i) > period/2
    phase_shift(i)=phase_shift(i)-period;
elseif phase_shift(i) < -period/2
    phase_shift(i)=phase_shift(i)+period;
end

   
end
Tp_new=0:time_interval:period*N_rep;
plot(Tp_new,phase_shift,'k*');
xlabel('Initial phase (hr)')
ylabel('Phase shifts (hr)')

function dadt=ode2(time,B,perturbation_time,P_light,N_light,P_time)

dadt=zeros(size(B));
M=B(1);Pt=B(2);

 if time>=perturbation_time & time<=perturbation_time+P_time
  a=P_light;
 else
    a=N_light;
 end

vm=1;km=0.1;vp=0.5;kp1=10;kp2=0.03;kp3=0.1;Pcrit=0.1;Jp=0.05;keq=200-a;

q=2/(1+sqrt(1+8*keq*Pt));

P1=q*Pt;
P2=keq*P1^2;

dadt(1)=(vm)/(1+(Pt*(1-q)/(2*(Pcrit)))^2)-km*M;
dadt(2)=vp*M-(kp1*Pt*q+kp2*Pt)/((Jp)+Pt)-kp3*Pt;
