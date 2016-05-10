% MAIN DRIVER
%% Initialize Parameters
% PK Parameters
mtfv = 287.2; %Molecular weight of tenofovir (g/mol)
V1= 244; %Volume of central compartment (L)
V2= 464.54; %Volume of peripheral compartment (L)
Vcell= 0.28*10^-12; %Volume of PBMC compartment
Fbio = 0.32; %Bioavailability of tenofovir
ka = 1; %Absorption constant (1/h)
Q = 71.41; %Intercompartment clearance (L/h)
Cl = 29.28; %Clearance (L/h)
kout = 0.0144375; %Elimination of TFV-DP from PMBC (1/h)
kcat1 = 2.4*3600; %Turnover number for AK2 and TFV (1/h)
kcat2 = 0.12*3600; %Turnover number for NDKA and TFV-MP
Km1 = 3*10^6; %Michaelis-Menten constant for AK2 and TFV (nM = nmol/L)
Km2 = 0.29*10^6; %Michaelis-Menten constant for NDKA and TFV-MP (nM)
Km = 770*10^3; %Michaelis-Menten constant for TFV uptake (nM)
% Km = 2*10^-4*Km; %Assumed for saturating conditions
Vmax = 1.77*10^9/mtfv; %Maximum rate of TFV uptake (nM/h)
E1 = 56; %Enzyme concentration for hAK2 (nM)
E2 = 287.7; %Enzyme concentration for NDA
GSHi = 3.1*10^6; %Intracellular GSH concentration (nM)
GSHe = 20*10^3; %Extracellular GSH concentration (nM)
fu = 0.93; %unbound fraction of TFV
Vbl = 5;
N = 0.9*10^6*10^3; %PBMC per ml of blood (cells/L)
% kleak = 9.2*10^3; % Leakage of TFV from literature
kleak = 3.4*10^10; % Leakage of TFV, assumed in the model
Vcell2 = Vcell*N*Vbl;
p = [0 V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu Vbl N kleak]';

% PD Parameters
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
IC_50 = 0.001*10^3; %nmol/L, range from literature is 0.04 - 8.5 umol/L
p_viral=[gammaT, gammaM, deltaT, deltaM, deltaPICT, deltaPICM, kT, kM, NhatT,...
    NhatM, NT,NM,deltaT1,deltaT2,deltaM1,deltaM2,CL_n,CL_in IC_50]';

% Initialize Tu,Mu,T1,M1,T2,M2,VI,VNI
y0_viral = set_initial_conditions(1500);

%% Part 2: Analyze PK/PD Model with Full Adherence
run('driver_part2.m');