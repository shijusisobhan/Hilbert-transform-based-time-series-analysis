# Hilbert transform based time series analysis of the circadian gene regulatory network
## Background
 Different tools are used to extract different properties of circadian rhythms. Instead of using different
techniques to obtain different properties of the circadian rhythms, this work propose the Hilbert transform (HT)-based numerical method which can be used to
extract most of the important circadian properties from the circadian time series. Here we providing Matlab code and data for generating figures appeared in the paper [Shiju, S., & Sriram, K. (2019). Hilbert transform based time series analysis of the circadian gene regulatory network, IET System Biology, ISSN : 1751-8849](https://digital-library.theiet.org/content/journals/10.1049/iet-syb.2018.5088)
## Tyson_ode.m

This is a function file that containe ODE of Tyson et.al (1999) model

This function is called by the programs named 'Arnold_tonguel.m','Arnold_onion.m' 'Devil_staircase.m', 
'period_dt_Tyson1999.m'.

This function requires 5 input arguments to execute and they are

t: Time

B:Initial points of the variable

a: amplitude of the forcing signal

tau: period of forcing signal

PP: photo period of the forcing signal

If it is free run (no external signal), then a=0;

## Period_det_Tyson1999.m

This is the program to determine the period of the tyson et al model using Hilbert transform. Using this program we can generate Fig 3.

##  Tyson_stochastic_period.m

This is the program to determine the period of the stochastic verion of tyson et al model using Hilbert transform.
Here we use a parameter 'Omg' for various number of molecules.
using this program we can generate Fig 4.

##  Period_sensitivity_HT.m

This is the program to determine the period sensitivity of the tyson model using HT. 
dP is the amount of change we made in the parameter (10% of of each parameter value).
for finding sensitivity of each parameter add dP to corresponding parameter value. eg. vm=1+dP
using this program we can generate Fig 5.

## PRC_Tyson_HT
This is the program to construc PRC of the tyson model using HT.
using this program we can generate Fig 8.

## PRC_Drosophila_goldbeter1999.m

This is the program to construc PRC of the Goldbeter Drosophila model using HT.
using this program we can generate Fig 7A.

## PRC_Goodwin_Rouff.m

This is the program to construc PRC of the Rouff Neurospora model using HT.
using this program we can generate Fig 7B.

## Devil_staircase.m
This is the program to construc PRC of the tyson model using HT.
using this program we can generate Fig 9A.

## Arnold_tongue.m
This is the program to construc Arnold's tongue of the tyson model using HT.
using this program we can generate Fig 9B.

## Tyson_phaseslip.m

This is the program to calculate the  phase slip of the tyson model using HT.
using this program we can generate Fig 11.

## Tyson_stochastic_phaseslip.m

This is the program to calculate the  phase slip of the stochastic simulation tyson model using HT.
using this program we can generate Fig 12.

## Slip_experimental.m

This is the program to calculate period and phase slip of the SCN experimental data. Data is 
provided in the same folder.


## Cite this work
If you are using these code in your work, you can cite as:
Shiju, S., & Sriram, K. (2019). Hilbert transform based time series analysis of the circadian gene regulatory network, IET System Biology, ISSN : 1751-8849
