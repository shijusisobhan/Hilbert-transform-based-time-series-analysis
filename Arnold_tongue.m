clc;clear all;warning off
fs=100;
t=0:1/fs:2500;

B=0.2*ones(1,2);% initial condition
a=0:0.1:0.6;%amplitude of forcing signal
tau=23:0.5:26;%period of forcing signal
PP=50;%photoperiod;


period1=[];
Entrainment_amplitude1=[];
period2=[];
Entrainment_amplitude2=[];
period3=[];
Entrainment_amplitude3=[];
h = waitbar(0,'Please wait...');
for i=1:length(tau);
    for j=1:length(a)
    [t,A]=ode45(@Tyson_ode,t,B,[],a(j),tau(i), PP);
    x=A(:,1);
    
x1=detrend(x);
z = hilbert(x1);
instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
period=1/mean(instfreq(100000:245000))


    
if abs(2*period-tau(i))<=0.1
    Entrainment_amplitude1(end+1)=a(j);
    period1(end+1)=tau(i);
    
    
    
elseif abs(period-tau(i))<=0.1
    Entrainment_amplitude2(end+1)=a(j);
    period2(end+1)=tau(i);
    
elseif abs(0.5*period-tau(i))<=0.1
    Entrainment_amplitude3(end+1)=a(j);
    period3(end+1)=tau(i);
end
    end
    waitbar(i/length(tau),h)
end
    plot(period1,Entrainment_amplitude1,'b.')
hold on
plot(period2,Entrainment_amplitude2,'r.')
hold on
plot(period3,Entrainment_amplitude3,'k.')