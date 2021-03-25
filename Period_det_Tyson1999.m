clc;clear all;
%% This is the program to calculate period of deterministic model (Tyson et al)
%% approximate time to run: 17sec
dt=0.01;
fs=1/dt;
t=0:dt:7000;

B=0.1*ones(1,2);% initial condition
a=0;%amplitude of forcing signal
tau=24;%period of forcing signal
PP=50;%photoperiod;
%% ************************************************
[t,A]=ode23s(@Tyson_ode,t,B,[],a,tau,PP);

%% period calculation using HT
x=A(:,1);
x1=detrend(x);
z = hilbert(x1);
phi=angle(z);
xp=phi;

zerIdx1=[];
for i=1:length(xp)-1
   if (abs(xp(i)-xp(i+1))>(1.5*pi))
      zerIdx1(end+1)=i; % save index of zero-crossing
   end
end

t1=t(zerIdx1); 
Period_switching=[];
for j=1:length(t1)-1
    period3(j)=t1(j+1)-t1(j);
    
    if period3(j)>5
       Period_switching(end+1)=period3(j);
    end
        
end

subplot(311)
plot(t,x)
axis([2000 2072 0 3])
xlabel('Time(h)')
ylabel('M')
subplot(312)
plot(t,phi);
axis([2000 2072 -4 4])
xlabel('Time(h)')
ylabel('\theta (rad)')
subplot(313)
hist(Period_switching(5:end-5),1)
axis([22 26 0 300])
xlabel('period(h)')
ylabel('occurrence')
