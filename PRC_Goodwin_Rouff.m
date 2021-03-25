function PRC_Goodwin_Rouff

%% Required data

time=0:0.01:4000;

number_of_points_required=24;%How many ponits required for PRC

N_rep=1;%Number of PRC cycle required

Number_cycle=10;%Number of cycle after which phase calculating

T_ref=4.3 ;%specify the circadian time of peak value of the fist varible

P_light=1;%Intensity of perturbation(for PRC)

N_light=0; %Unpertubed intensity

P_time=1;%perturbation duration

B=[0.039 0.18 1.7];
%% 
option = odeset('RelTol', 1e-8);

%[t1,A1]=ode15s(@ode1,t,initial_condition,option,perturbation time);

[t1,A1]=ode23s(@ode2,time,B,option,[],N_light,N_light,[]);

%******************************************************************************
%% Period claculation

x1=find(time>3000);
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

kk1=find(time>=3000 & time<=3000+period);

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
[time,A]=ode23s(@ode2,time,B,option,perturbation_time,P_light,N_light,P_time);

T_phase_calculation=T_initial+Number_cycle*period;%Phase calculation after N cycle

%T1=find(time>=T_phase_calculation & time<=T_phase_calculation+period*N_rep);%phase shift calculate after 20 cycle
T1=find(time>=T_phase_calculation);
T_new=time(T1);

X=A1(:,1);Y=A(:,1);
x_pt=X(T1);y_pt=Y(T1);

%%  phase difference calculation using hilbert transform

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

T_N=data(:,5);
P_wild=data(:,6);

%%% Plot PRC:

plot(T_N,P_wild,'ro')

legend('Simulation','experiment')


function dadt=ode2(time,B,perturbation_time,P_light,N_light,P_time)

dadt=zeros(size(B));

X=B(1);Y=B(2);Z=B(3);

k1=1;k2=1;k3=1;k4=0.2;k5=0.2;k6=0.1;

if time>=perturbation_time & time<=perturbation_time+P_time
  light=P_light;
 else
    light=N_light;
 end

dadt(1)=(light+k1)/(1+Z^9)-k4*X;
dadt(2)=k2*X-k5*Y;
dadt(3)=k3*Y-k6*Z;