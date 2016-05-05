function out = virus_dynamics(etaterm)

gammaT=2e9; % birth rate of uninfected T-cells, 1/day
gammaM=6.9e7; % birth rate of uninfected macrophages, 1/dzy
deltaT=0.02; % death rate of uninfected T-cells, 1/day
deltaM=0.0069; % death rate of uninfected macrophages

% intracellular degradation of essential components of the pre-integration 
% complex, e.g., by the host cell proteasome, which return early infected 
% T-cells and macrophages to an uninfected stage, respectively
deltaPICT=0.35;
deltaPICM=0.0035;

% rate constants of proviral integration into the host cell’s genome
kT=0.35;
kM=0.07; 

% total number of released infectious and non-infectious 
NhatT=1000;
NhatM=100;
% virus from late infected T-cells and macrophages
% rates of release of infectious virus
NT=0.67*NhatT;
NM=0.67*NhatM;
% death rate constants of T1, T2, M1,M2
deltaT1=deltaT;
deltaT2=1;
deltaM1=deltaM;
deltaM2=0.09; 
CL_n=2.3; % clearance rate of free virus by the immune system
CL_in=23;

y0 = set_initial_conditions(24);
p=[gammaT, gammaM, deltaT, deltaM, deltaPICT, deltaPICM, kT, kM, NhatT,...
    NhatM, NT,NM,deltaT1,deltaT2,deltaM1,deltaM2,CL_n,CL_in,etaterm];
TimeLen = 24; % hours
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode15s(@virus_dynamics_eqns,[0 TimeLen],y0,options,p);
end



