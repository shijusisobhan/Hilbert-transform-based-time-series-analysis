clc;clear all;
warning off;

%% This is the program for estimating the period of stochstic simulation for tyson et al. model

%% initialization
%%%general simulation parameters
t=0; %start time 
t_end=2500; %end time
t_sample=1; %sample interval for gathering data
k=1; %counter for waitbar update 
alpha=10^2; %parameter for updating waitbar (increase for shorter runtime)

%%%model parameters
%initial numbers for each chemical species 
M = 1; %mRNA
Pt = 5; %protien

vm=1;km=0.1;vp=0.5;kp1=10;kp2=0.03;kp3=0.1;Pcrit=0.1;Jp=0.05;keq=200;

 
% keq=200;

%%%arrays to store results
j=1; %counter for arrays
t_array(1,t_end/t_sample+1)=0; t_array(1,j)=t; %time array and initial value
M_array(1,t_end/t_sample+1)=0; M_array(1,j)=M; %M array and initial value
Pt_array(1,t_end/t_sample+1)=0; Pt_array(1,j)=Pt; %Pt array and initial value


Omg=1e2;% Number of molecules

%% SSA
tic %start timing the Gillespie loop
w=waitbar(0,'running SSA...');
while t < t_end,
    %calculate rxn propensities       
        h=[Omg*vm/(1+(Pt*(1-(2/(1+sqrt(1+8*(keq/Omg)*Pt))))/(2*Omg*Pcrit))^2) km*M...
            vp*M (Omg*kp1*Pt*(2/(1+sqrt(1+8*(keq/Omg)*Pt)))+Omg*kp2*Pt)/(Jp*Omg+Pt) kp3*Pt];      
    %combined rxn hazard
    h0 = sum(h);   
    %calculate time to next event
    r1=rand;
    while r1 == 0,
        r1=rand;
    end
    t_next = ((1/h0)*(log(1/r1)));
    
    %update time
    t = t + t_next;   
    %determine next reaction
    i=1; mu=0; amu=0; r2=rand;
    while amu < r2*h0,
        mu = mu + 1;
        amu = amu + h(i); 
        i = i + 1;
    end
    
    %reactions
    if mu == 1 %Reaction 1
        M = M + 1;
    elseif mu == 2 %Reaction 2
        M = M - 1;
    elseif mu == 3 %Reaction 3
        Pt = Pt + 1;
    elseif mu == 4 %Reaction 4
        Pt = Pt - 1;
        elseif mu == 5 %Reaction 5
        Pt = Pt - 1;
          end
    
    %store/output time and species
    if t >= j*t_sample
        j=j+1;
        t_array(1,j)=j;
        M_array(1,j)=M;
        Pt_array(1,j)=Pt;
       
    end    
    
    %update waitbar
    if t >= k*alpha*t_sample
        k=k+1;
        waitbar(t/t_end)
    end
end 
close(w)
toc


x2=M_array;
t=t_array;
x1=detrend(x2);  
z = hilbert(x1);
Fs=1;
fs=1;
phi=(angle(z));

%period calculation method1: positive mean crossing of M
xp=M_array;
N=mean(xp);
zerIdx=[];
 
for i=1:length(xp)-1
   if ((xp(i)<N && xp(i+1)>N))
      zerIdx(end+1)=i; % save index of zero-crossing
   end
end
t1=t(zerIdx); 

Period_zero_crossing=[];
for j=1:length(t1)-1
    period2(j)=t1(j+1)-t1(j);
    
%     if period2(j)>5
       Period_zero_crossing(end+1)=period2(j);
%     end
end
%  hist(Period_zero_crossing(5:end-5))


%% period calculation method2: switching time of phase (HT)
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
       Period_switching(end+1)=period3(j);    
end

%% **************************************************************
subplot(411)
plot(t_array,x2)
xlabel('time(h)')
ylabel('M')
subplot(412)
plot(t_array,phi)
xlabel('time(h)')
ylabel('\theta(rad)')
subplot(413)
hist(Period_zero_crossing(5:end-5))
title('mean crossing method')
xlabel('period(h)')
ylabel('occurence')
subplot(414)
hist(Period_switching(5:end-5))
xlabel('period(h)')
ylabel('occurence')
title('HT')
% M1=mean(Period_zero_crossing(5:end-5))
% SD1=std(Period_zero_crossing(5:end-5))
% M2=mean(Period_switching(5:end-5))
% SD2=std(Period_switching(5:end-5))