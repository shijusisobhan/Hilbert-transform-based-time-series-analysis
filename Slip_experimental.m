clc;clear all;warning off;

%% program for calculating experimental data

%% Experimetal data collected from 
%Data Source: Bordyugov, Grigory, et al. "Tuning the phase of circadian entrainment." 
%Journal of The Royal Society Interface 12.108 (2015): 20150282.							



%% For T24 temperature cycle

A=xlsread('SCN_Entrainment_T24.xlsx');
T_E=24;t=A(:,1);xs=A(:,2);LD=A(:,3);

%% For T26 temperature cycle

% A=xlsread('SCN_Entrainment_T26_T22.xlsx');
% T_E=26;t=A(:,1);xs=A(:,2);LD=A(:,3);

%% For T22 temperature cycle

%  A=xlsread('SCN_Entrainment_T26_T22.xlsx');
%   T_E=22;t=A(:,1);xs=A(:,4);LD=A(:,5);

%% ****************************************************
dt=1/12;
t_end=t(end);
%% *******************************
x2=xs;
x1=detrend(x2);  
z = hilbert(x1);
phi=(angle(z));
% period calculation
x=phi;
N=0;
zerIdx=[];
for i=1:length(x)-1
   if ((x(i)>N && x(i+1)<N))
      zerIdx(end+1)=i; % save index of phase switching
   end
end
t1=t(zerIdx); 
for j=1:length(t1)-1
    period1(j)=t1(j+1)-t1(j);
end
period1;
%% period calculation (cycle by cycle)
period_bin=period1(1)*ones(1,(t1(1)/dt));% initilizing the Instantaneous period
% 
for td=1:length(period1)
    ynm=period1(td)*ones(1,round((period1(td)/dt)));
   period_bin=[period_bin ynm];   
end

period_bin=period_bin';
slip(1)=0;

for i=1:length(period_bin);
    slip(i+1)=slip(i)+dt*((period_bin(i)-T_E)./period_bin(i));
end


%% Ploting commands

subplot(411)
plotyy(t/24,xs,t/24,LD);
grid on
ylabel('Bioluminescence')

subplot(412)
plot(t/24,phi)
grid on
ylabel('\theta');

subplot(413)
plot((t(1:length(period_bin)))/24,period_bin)
grid on
ylabel('T_{inst}');

subplot(414)
 plot(t(1:length(slip))/24,slip/24);
 grid on
 ylabel('\psi');