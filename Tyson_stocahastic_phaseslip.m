clc;clear all;
warning off;
%% initialization
%%%general simulation parameters
t=0; %start time 
t_end=8000; %end time
t_sample=1; %sample interval for gathering data
k=1; %counter for waitbar update 
alpha=10^2; %parameter for updating waitbar (increase for shorter runtime)

%%%model parameters
%initial numbers for each chemical species 
M = 1; %MRNA
Pt = 5; %protien

vm=1;km=0.1;vp=0.5;kp1=10;kp2=0.03;kp3=0.1;Pcrit=0.1;Jp=0.05;

 
% keq=200;

%%%arrays to store results
j=1; %counter for arrays
t_array(1,t_end/t_sample+1)=0; t_array(1,j)=t; %time array and initial value
M_array(1,t_end/t_sample+1)=0; M_array(1,j)=M; %promoter array and initial value
Pt_array(1,t_end/t_sample+1)=0; Pt_array(1,j)=Pt; %mRNA array and initial value


Omg=1000;
 a=0.15;
 T_E=23;
%% SSA
tic %start timing the Gillespie loop
w=waitbar(0,'running SSA...');
while t < t_end,
 keq=(-100*a*square(2*pi*(1/T_E)*t,50)+100*a)+(200-200*a);
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
    if mu == 1 %transcription
        M = M + 1;
    elseif mu == 2 %mRNA decay
        M = M - 1;
    elseif mu == 3 %translation 
        Pt = Pt + 1;
    elseif mu == 4 %protein decay
        Pt = Pt - 1;
        elseif mu == 5 %protein decay
        Pt = Pt - 1;
          end
    
    %store/output time and species
    if t >= j*t_sample
        j=j+1;
        t_array(1,j)=j;
        M_array(1,j)=M;
        Pt_array(1,j)=Pt;
        t_steps(1,j)=t_next;
       
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

x1=detrend(x2);  
z = hilbert(x1);
Fs=1;
fs=1;
phi=(angle(z));

% instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
%  period1=1/mean(instfreq(1000:3500));

% period calculation
x=phi;
N=0;
zerIdx=[];
for i=1:length(x)-1
   if ((x(i)>N && x(i+1)<N))
      zerIdx(end+1)=i; % save index of zero-crossing
   end
end
t1=t_array(zerIdx); 
for j=1:length(t1)-1
    period1(j)=t1(j+1)-t1(j);
end
period1;
M1=mean(period1);%find out the mean period
Period_of_oscillation=M1

% %% period calculation (cycle by cycle)
cycle=round(t_end/T_E);
period_bin=[];
% 
for td=1:cycle-20
    ynm=period1(td)*ones(1,(T_E/1));
  
   period_bin=[period_bin ynm];   
end

step_size=1+t_steps';
slip(1)=0;
for i=1:length(period_bin);
    slip(i+1)=slip(i)+step_size(i)*((period_bin(i)-T_E)./period_bin(i));
end

 plot(t_array(1:length(slip)),slip/24);

