%% Driver
clear all;
close all;

%% Initialize Values
mtfv = 287.2*10^3; %Molecular weight of tenofovir (mg/mol)

dose = 300; %300 mg, taken orally
V1= 244; %Volume of central compartment (L)
V2= 464.54; %Volume of peripheral compartment (L)
Vcell= 0.28*10^-12; %Volume of PBMC compartment
Fbio = 0.32; %Bioavailability of tenofovir
ka = 1; %Absorption constant (1/h)
Q = 71.41; %Intercompartment clearance (L/h)
Cl = 29.28; %Clearance (L/h)
kout = 0.0144; %Elimination of TFV-DP from PMBC (1/h)
kcat1 = 2.4*3600; %Turnover number for AK2 and TFV (1/h)
kcat2 = 0.12*3600; %Turnover number for NDKA and TFV-MP
Km1 = 3*10^-3; %Michaelis-Menten constant for AK2 and TFV (M)
Km2 = 0.29*10^-3; %Michaelis-Menten constant for NDKA and TFV-MP (M)
Km = 770*10^-6; %Michaelis-Menten constant for TFV uptake (M)
Vmax = 1.77/mtfv; %Maximum rate of TFV uptake (M/h)
E1 = 56*10^-9; %Enzyme concentration for hAK2 (M)
E2 = 287.7*10^-9; %Enzyme concentration for NDA
GSHi = 3.1*10^-3; %Intracellular GSH concentration
GSHe = 20*10^-6; %Extracellular GSH concentration
fu = 0.93; %unbound fraction of TFV
F = 4.178*10^-5; %unionized fraction of TFV
Vbl = 5;
N = 0.9*10^6*10^3; %PBMC per ml of blood (cells/L)
kleak = 9.2*10^3; % Leakage of TFV

p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';

%% Run Simulation
TimeLen = 24;
OutputVar = 1:3;
[AUC,Balance,t,y] = Tenofovir(p,OutputVar,TimeLen);
figure;
ax1=subplot(2,2,1);
plot(ax1,t,y(:,1),'k',t,y(:,2),'r.',t,y(:,3),'b-','linewidth',3)
title(ax1,'Concentration of Drug in Compartments')
ylabel(ax1,'[D] (M)')
xlabel(ax1,'time (hrs)')

ax4=subplot(2,2,4);
plot(ax4,T1,Balance,'k','linewidth',3)
title(ax4,'Molecular Balance') %(zero = balance)
ylabel(ax4,'Balance of Drug (mol)')
xlabel(ax4,'time (hrs)')