%% Driver
clear all;
close all;


%% Initialize PK parameters
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
F = 4.178*10^-5; %unionized fraction of TFV
Vbl = 5;
N = 0.9*10^6*10^3; %PBMC per ml of blood (cells/L)
% kleak = 9.2*10^3; % Leakage of TFV
kleak = 3.4*10^10; % Leakage of TFV
Vcell2 = Vcell*N*Vbl;
p = [0 V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';

%% Initialize PD parameters: Viral Dynamics
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

p_viral=[gammaT, gammaM, deltaT, deltaM, deltaPICT, deltaPICM, kT, kM, NhatT,...
    NhatM, NT,NM,deltaT1,deltaT2,deltaM1,deltaM2,CL_n,CL_in]';
%% Run Simulation
dose_set = [75, 150, 300, 600];
TimeLen = 15*24;
OutputVar = 1:15;
y0_viral = set_initial_conditions(400);
for i = 1:4
    dose = dose_set(i)*10^6; %dose in ng, taken orally
    p(1) = dose;
    [AUC,Balance,t,y] = Tenofovir(p,p_viral,y0_viral,OutputVar,TimeLen);
    yset{i} = y;
    tset{i} = t;
    balanceset{i} = Balance;
end
y1 = yset{1};
y2 = yset{2};
y3 = yset{3};
y4 = yset{4};
%% Plot Figures
figure;
ax1=subplot(2,2,1);
plot(ax1,tset{1},y1(:,1)/(V1*10^3),'k',tset{2},y2(:,1)/(V1*10^3),'r',tset{3},y3(:,1)/(V1*10^3),'b',tset{4},y4(:,1)/(V1*10^3),'g','linewidth',3)
title(ax1,'Concentration of TFV in Central Compartment')
ylabel(ax1,'TFV (nmol/mL)')
xlabel(ax1,'time (hrs)');
ax4=subplot(2,2,2);
plot(ax4,tset{1},y1(:,4)/(Vcell2*10^3),'k',tset{2},y2(:,4)/(Vcell2*10^3),'r',tset{3},y3(:,4)/(Vcell2*10^3),'b',tset{4},y4(:,4)/(Vcell2*10^3),'g','linewidth',3)
title(ax4,'Concentration of TFV-MP in PBMC') %(zero = balance)
ylabel(ax4,'TFV-MP (nmol/mL)')
xlabel(ax4,'time (hrs)')

ax4=subplot(2,2,3);
plot(ax4,tset{1},y1(:,5)/(Vcell2*10^3),'k',tset{2},y2(:,5)/(Vcell2*10^3),'r',tset{3},y3(:,5)/(Vcell2*10^3),'b',tset{4},y4(:,5)/(Vcell2*10^3),'g','linewidth',3)
title(ax4,'Concentration of TFV-DP in PBMC') %(zero = balance)
ylabel(ax4,'TFV-DP (nmol/mL)')
xlabel(ax4,'time (hrs)')

ax4=subplot(2,2,4);
plot(ax4,tset{1},balanceset{1},'k',tset{2},balanceset{2},'r',tset{3},balanceset{3},'b',tset{4},balanceset{4},'g','linewidth',3)
title(ax4,'Molecular Balance') %(zero = balance)
ylabel(ax4,'Balance of Drug (nmol)')
xlabel(ax4,'time (hrs)')
legend('75 mg','150 mg','300 mg','600 mg');
% figure;
% ax1=subplot(1,3,1);
% plot(ax1,t,log10(y(:,12)),'linewidth',3)
% title(ax1,'Number of infected T-cells after proviral genomic integration')
% ylabel(ax1,'log10 T2')
% xlabel(ax1,'time (hrs)')
% 
% 
% ax4=subplot(1,3,2);
% plot(ax4,t,log10(y(:,13)),'linewidth',3)
% title(ax4,'Number of infected macrophages after proviral genomic integration')
% ylabel(ax4,'log10 M2')
% xlabel(ax4,'time (hrs)')
% 
% ax4=subplot(1,3,3);
% plot(ax4,t,log10(-(y(:,14)+y(:,15) - (y0_viral(7)+y0_viral(8)))),'linewidth',3)
% title(ax4,'Viral Load Decay')
% ylabel(ax4,'Log10 Viral Load Decay')
% xlabel(ax4,'time (hrs)')
