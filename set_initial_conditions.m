% compute and return initial conditions of the system:
function y0 = set_initial_conditions(T)
gammaT=(2*10^9)/24; % birth rate of uninfected T-cells,1/day to 1/hour
gammaM=(6.9*10^7)/24; % birth rate of uninfected macrophages, 1/day to 1/hour
deltaT=0.02/24; % death rate of uninfected T-cells and T1, 1/day to 1/hour
deltaM=0.0069/24; % death rate of uninfected macrophages and M1

% intracellular degradation of essential components of the pre-integration
% complex, e.g., by the host cell proteasome, which return early infected
% T-cells and macrophages to an uninfected stage, respectively
deltaPICT=0.35/24;
deltaPICM=0.0035/24;

% rate constants of proviral integration into the host cell’s genome
kT=0.35/24;
kM=0.07/24;

% total number of released infectious and non-infectious
NhatT=1000/24;
NhatM=100/24;
% virus from late infected T-cells and macrophages
% rates of release of infectious virus
NT=0.67*NhatT;
NM=0.67*NhatM;
% death rate constants of T1, T2, M1,M2
deltaT1=0.5/24;
deltaT2=1.4/24;
deltaM1=deltaM;
deltaM2=0.09/24;
CL_n=2.3/24; % clearance rate of free virus in uninfected individuals
CL_in=23/24; % clearance rate of free virus in infected individuals

p=[gammaT, gammaM, deltaT, deltaM, deltaPICT, deltaPICM, kT, kM, NhatT,...
    NhatM, NT,NM,deltaT1,deltaT2,deltaM1,deltaM2,CL_n,CL_in];
% pseudo-initial loads (obtained by pure deterministic simulations, used as initial conditions for pre-run)
Initial_T  = 0.0361e10;
Initial_T1 = 0.5386e10;
Initial_T2 = 1.885e10;
Initial_M  = 0.3485e10;
Initial_M1 = 0.0581e10;
Initial_M2 = 0.0451e10;
Initial_VI  = 5.5885e10;
Initial_VNI = 2.7694e10;
% set pseudo-initial species numbers
y0 = [Initial_T;Initial_M;Initial_T1;Initial_M1;Initial_T2;Initial_M2;Initial_VI;Initial_VNI];

% compute initial species numbers
TimeLen = T; % hours
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[t,y] = ode15s(@virus_dynamics_eqns,[0 TimeLen],y0,options,p);

% set initial species vector
y0 = round(y(end,:));
% 
% figure;
% ax1=subplot(1,3,1);
% plot(ax1,t,y(:,5),'linewidth',3)
% title(ax1,'Number of infected T-cells after proviral genomic integration')
% ylabel(ax1,'T2')
% xlabel(ax1,'time (hrs)')
% 
% 
% ax4=subplot(1,3,2);
% plot(ax4,t,y(:,6),'linewidth',3)
% title(ax4,'Number of infected macrophages after proviral genomic integration')
% ylabel(ax4,'M2')
% xlabel(ax4,'time (hrs)')
% 
% ax4=subplot(1,3,3);
% plot(ax4,t,y(:,7)+y(:,8),'linewidth',3)
% title(ax4,'Viral Load') %(zero = balance)
% ylabel(ax4,'VI + VNI')
% xlabel(ax4,'time (hrs)')
end
