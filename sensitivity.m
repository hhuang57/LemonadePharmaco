%% Initialize PK parameters
mtfv = 287.2; %Molecular weight of tenofovir (g/mol)
OutputVar = 1:5;
TimeLen=120;
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
dose=300*10^6; %dose in ng, taken orally
%%
TimeLen=1/3;
ParamDelta = 0.05;
p0= [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';
OutputVar = 1:5;
[auc,Balance0,t0,y0] = Tenofovirpop(p0,OutputVar,TimeLen);
auc0=auc;
for i=1:length(p0)
    p=p0;
    p(i)=p0(i)*(1.0+ParamDelta);
    [auc,~,t1,y1] = Tenofovirpop(p,OutputVar,TimeLen);
    auc(i)=auc;
    Sens(i) = ((auc(i)-auc0)/auc0)/((p(i)-p0(i))/p0(i));
    SensAbs(i) = ((auc(i)-auc0)/auc0);
end

figure;
bar(Sens);
ylim([-1,1]);
title('Local Sensitivity: AUC of TFV');
