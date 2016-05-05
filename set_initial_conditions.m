% compute and return initial conditions of the system:
function y0 = set_initial_conditions(T)
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

p=[gammaT, gammaM, deltaT, deltaM, deltaPICT, deltaPICM, kT, kM, NhatT,...
    NhatM, NT,NM,deltaT1,deltaT2,deltaM1,deltaM2,CL_n,CL_in];
% pseudo-initial loads (obtained by pure deterministic simulations, used as initial conditions for pre-run)
Initial_T  = 2.0e11 + 0.02*2.0e11 - 1.0e6 - 1.0e5;
Initial_T1 = 1.0e5;
Initial_T2 = 1.0e6;
Initial_M  = 2.0e9 + 0.02*2.0e9 - 1.0e6 - 1.0e5;
Initial_M1 = 1.0e3;
Initial_M2 = 1.0e4;
Initial_VI  = 17.5e6;
Initial_VNI = 5e6;
% set pseudo-initial species numbers
y0 = [Initial_T;Initial_M;Initial_T1;Initial_M1;Initial_T2;Initial_M2;Initial_VI;Initial_VNI];

% compute initial species numbers
TimeLen = T; % hours
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[t,y] = ode15s(@virus_dynamics_eqns,[0 TimeLen],y0,options,p);

% set initial species vector
y0 = round(y(end,:));

figure;
ax1=subplot(1,3,1);
plot(ax1,t,y(:,3),'linewidth',3)
title(ax1,'Number of infected T-cells after proviral genomic integration')
ylabel(ax1,'T2')
xlabel(ax1,'time (hrs)')


ax4=subplot(1,3,2);
plot(ax4,t,y(:,4),'linewidth',3)
title(ax4,'Number of infected macrophages after proviral genomic integration')
ylabel(ax4,'M2')
xlabel(ax4,'time (hrs)')

ax4=subplot(1,3,3);
plot(ax4,t,y(:,7)+y(:,8),'linewidth',3)
title(ax4,'Viral Load') %(zero = balance)
ylabel(ax4,'VI + VNI')
xlabel(ax4,'time (hrs)')
end
