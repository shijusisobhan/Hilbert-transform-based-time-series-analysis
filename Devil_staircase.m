clc;clear all;warning off
fs=100;
t=0:1/fs:2500;

B=0.2*ones(1,2);% initial condition
a=0.5;%amplitude of forcing signal
tau=0:0.5:60;%period of forcing signal
PP=50;%photoperiod;

for i=1:length(tau);
    [t,A]=ode23s(@Tyson_ode,t,B,[],a,tau(i),PP);
     x=A(:,1);
     tau(i)

x1=detrend(x);
z = hilbert(x1);
instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
period1=1/mean(instfreq(100000:245000))
%% **************************************************************
period_number(i)=period1/tau(i);
end
period_number';

plot(1./tau,period_number);
xlabel('forcing frequency')
ylabel('Period number')