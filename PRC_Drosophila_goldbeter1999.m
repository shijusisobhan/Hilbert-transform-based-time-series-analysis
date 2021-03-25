function PRC_Drosophila_goldbeter1999

%% Required data

time=0:0.01:1500;

number_of_points_required=24;%How many ponits required for PRC

N_rep=1;%Number of PRC cycle required

Number_cycle=10;%Number of cycle after which phase calculating

T_ref=6;%specify the circadian time of peak value of the fist varible

%P_light=2;%Intensity of perturbation(for PRC)

P_light=1; %Intensity of perturbation (for VRC)

N_light=0; %Unpertubed intensity

P_time=3;%perturbation duration

B=0.1*ones(1,10);%Initial condition
%% 

option = odeset('RelTol', 1e-8);

%[t1,A1]=ode15s(@ode1,t,initial_condition,option,perturbation time);

[t1,A1]=ode45(@ode2,time,B,option,[],N_light,N_light,[]);

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

 h = waitbar(0,'Please wait...'); 
for i=1:length(Tp)
  perturbation_time=Tp(i);
[time,A]=ode45(@ode2,time,B,option,perturbation_time,P_light,N_light,P_time);

T_phase_calculation=T_initial+Number_cycle*period;%Phase calculation after N cycle

%T1=find(time>=T_phase_calculation & time<=T_phase_calculation+period*N_rep);%phase shift calculate after 20 cycle
T1=find(time>=T_phase_calculation);
T_new=time(T1);

X=A1(:,1);Y=A(:,1);
x_pt=X(T1);y_pt=Y(T1);

%%  phase difference calculation using hilbert transform

% 
h1=hilbert(detrend(x_pt));
h2=hilbert(detrend(y_pt));

p1 = unwrap(angle(h1)); p2 =unwrap(angle(h2)); % Instantaneous phase (wrapped)
tau=period/(2*pi);

p11=(p2-p1);

p=mean(p11(5*(24/0.01):end-2*(24/0.01)));

 phase_shift(i)=p*tau;
 
 if phase_shift(i) > period/2
    phase_shift(i)=phase_shift(i)-period;
elseif phase_shift(i) < -period/2
    phase_shift(i)=phase_shift(i)+period;
end
  

 waitbar(i/length(Tp),h)
end

Tp_new=0:time_interval:period*N_rep;

plot(Tp_new*(24/period),phase_shift,'b');
xlabel('Initial phase (hr)')
ylabel('Phase shifts (hr)')

hold on
data=xlsread('PRC_EXP.xlsx');%for plotting Experimental data 

T_N=data(:,1);
P_wild=data(:,2);

%%% Plot PRC:

plot(T_N,P_wild,'ro')

legend('Simulation','experiment')

function dadt=ode2(time,B,perturbation_time,P_light,N_light,P_time)


dadt=zeros(size(B));
Mp=B(1);P0=B(2);P1=B(3);P2=B(4);Mt=B(5);T0=B(6);T1=B(7);T2=B(8);C=B(9);Cn=B(10);

vsT=1;VmP=0.7;vmT=0.7;KmP=0.2;KmT=0.2;ksP=0.9;ksT=0.9;vdP=2;vdT=2;k1=0.6;k2=0.2;k3=1.2;
k4=0.6;KIP=1;KIT=1;KdP=0.2;KdT=0.2;n=4;K1P=2;K1T=2;K2P=2;K2T=2;K3P=2;K3T=2;K4P=2;K4T=2;kd=0.01;kdC=0.01;
kdN=0.01;V1P=8;V1T=8;V2P=1;V2T=1;
V3P=8;V3T=8;V4P=1;V4T=1;

VsP=1;

if time>=perturbation_time & time<=perturbation_time+P_time
  light=P_light;
 else
    light=N_light;
 end


dadt(1)=VsP*KIP^n/(KIP^n+Cn^n)-VmP*Mp/(KmP+Mp)-kd*Mp;

dadt(2)=ksP*Mp-V1P*P0/(K1P+P0)+V2P*P1/(K2P+P1)-kd*P0;

dadt(3)=V1P*P0/(K1P+P0)-V2P*P1/(K2P+P1)-V3P*P1/(K3P+P1)+V4P*P2/(K4P+P2)-kd*P1;

dadt(4)=V3P*P1/(K3P+P1)-V4P*P2/(K4P+P2)-k3*P2*T2+k4*C-vdP*P2/(KdP+P2)-kd*P2;

dadt(5)=vsT*KIT^n/(KIT^n+Cn^n)-vmT*Mt/(KmT+Mt)-kd*Mt;

dadt(6)=ksT*Mt-V1T*T0/(K1T+T0)+V2T*T1/(K2T+T1)-kd*T0;

dadt(7)=V1T*T0/(K1T+T0)-V2T*T1/(K2T+T1)-V3T*T1/(K3T+T1)+V4T*T2/(K4T+T2)-kd*T1;

dadt(8)=V3T*T1/(K3T+T1)-V4T*T2/(K4T+T2)-k3*P2*T2+k4*C-(vdT+light)*T2/(KdT+T2)-kd*T2;

dadt(9)=k3*P2*T2-k4*C-k1*C+k2*Cn-kdC*C;

dadt(10)=k1*C-k2*Cn-kdN*Cn;