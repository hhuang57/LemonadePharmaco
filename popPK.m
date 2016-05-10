%% Initialize PK parameters
mtfv = 287.2; %Molecular weight of tenofovir (g/mol)
OutputVar = 1:5;
TimeLen=240;
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
%p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';

%% Population

%things for ka
%http://www.cdc.gov/hiv/statistics/overview/
agedistrib=[0.1054,0.1957,0.2269,0.383,0.089]; %15-19,25-29,35-39,45-49,55-59
cumagedistrib=[0.089,0.383+0.089,0.2269+0.383+0.089,...
    0.1957+0.2269+0.383+0.089,0.1054+0.1957+0.2269+0.383+0.089]; 

%things for CL
BWstdev = [25.46,40.82,40.82,29.45,29.45]*0.454; %kg
BWmean = [96.8,169.6,169.6,173.2,173.2]*0.454; %kg
scrstdev=19; %muM
scrmean=83; %muM
number=36.2;%29.28/Fbio-136*BWmean/scrmean;

%things for kout
lb=.002;
ub=.026;
koutmean=(ub+lb)/2;
koutstdev=(ub-lb)/6; % ub and lb both 3 stdev away from the mean

%another thing for Cl
%Thirty-one percent of patients received lopinavir/ritonavir (LPV/r) as 
%a cotreatment, the addition of LPV/r treatment decreased tenofovir clearance by 25%
%Bouazza, Naïm, et al. "Population Pharmacokinetics of Tenofovir in HIV-1–Infected Pediatric Patients." JAIDS Journal of Acquired Immune Deficiency Syndromes 58.3 (2011): 283-288.

npeople=50;
info=zeros(4,npeople); %age,weight,scr,kout

for p1=1:npeople
    a=rand();
    if a<cumagedistrib(1)
        info(1,p1)=17;
    elseif a<cumagedistrib(2)
        info(1,p1)=27;
    elseif a<cumagedistrib(3)
        info(1,p1)=37;
    elseif a<cumagedistrib(4)
        info(1,p1)=47;
    else
        info(1,p1)=57;
    end
end

for p2=1:npeople
    group=(info(1,p2)-7)/10;
    info(2,p2)= BWstdev(group).*randn(1,1) + BWmean(group);
end

for p3=1:npeople
   info(3,p3)= scrstdev.*randn(1,1)+scrmean;
end

for p4=1:npeople
    info(4,p4)=koutstdev.*randn(1,1)+koutmean;
end

figure;
axx1=subplot(2,2,1);
hist(info(1,:));
title('Age Distribution');
axx2=subplot(2,2,2);
hist(info(2,:));
title('Weight Distribution');
axx3=subplot(2,2,3);
hist(info(3,:));
title('SCR Distribution');
axx4=subplot(2,2,4);
hist(info(4,:));
title('kout Distribution');
%% popPK for this population
aucall=zeros(4,npeople);
for sc=1:4
    for dum=1:npeople
        switch sc
            case 1
                ka= exp(-info(1,npeople)/50);
                Cl=29.28;
                kout = 0.0144375;
                p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';
                [outAUC,outBalance,outT,outY] = Tenofovirpop(p,OutputVar,TimeLen); 
                aucall(1,dum)=outAUC; 
            case 2
                ka=1;
                Cl=Fbio*(number+136*info(2,dum)/info(3,dum));
                kout = 0.0144375;
                p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';
                [outAUC,outBalance,outT,outY] = Tenofovirpop(p,OutputVar,TimeLen); 
                aucall(2,dum)=outAUC; 
            case 3
                ka=1;
                Cl=29.28;
                kout=info(4,dum);
                p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';
                [outAUC,outBalance,outT,outY] = Tenofovirpop(p,OutputVar,TimeLen); 
                aucall(3,dum)=outAUC; 
            case 4
                ka= exp(-info(1,npeople)/120);
                Cl=Fbio*(number+136*info(2,dum)/info(3,dum));
                kout=info(4,dum);
                p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';
                [outAUC,outBalance,outT,outY] = Tenofovirpop(p,OutputVar,TimeLen); 
                aucall(4,dum)=outAUC; 
        end
    end
end
figure;
boxplot(aucall');
ylabel('AUC of TFV-DP (hr-nmol/L)');
% title('Varying ka');
% xlabel('AUC of TFV (hr-nmol/L)');
%% cmax cmin
cminall=zeros(4,npeople);
cmaxall=zeros(4,npeople);
for sc=1:4
    for dum=1:npeople
        switch sc
            case 1
                ka= exp(-info(1,npeople)/50);
                Cl=29.28;
                kout = 0.0144375;
                p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';
                [outAUC,outBalance,outT,outY] = Tenofovir_multipledose(p,OutputVar,TimeLen); 
%                 plot(outT,outY);
%                 title('Dosing Every 24 Hours for 10 Days: Case 1');
%                 xlabel('Time (hr)');
%                 ylabel('TFV in Central Compartment (nmol/L)');
                cmin=min(outY(end-100:end));
                cmax=max(outY(end-100:end));
                cminall(sc,dum)=cmin;
                cmaxall(sc,dum)=cmax;
            case 2
                ka=1;
                Cl=Fbio*(number+136*info(2,dum)/info(3,dum));
                kout = 0.0144375;
                p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';
                [outAUC,outBalance,outT,outY] = Tenofovir_multipledose(p,OutputVar,TimeLen); 
%                 plot(outT,outY);
%                 title('Dosing Every 24 Hours for 10 Days: Case 2');
%                 xlabel('Time (hr)');
%                 ylabel('TFV in Central Compartment (nmol/L)');
                cmin=min(outY(end-100:end));
                cmax=max(outY(end-100:end));
                cminall(sc,dum)=cmin;
                cmaxall(sc,dum)=cmax;
            case 3
                ka=1;
                Cl=29.28;
                kout=info(4,dum);
                p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';
                [outAUC,outBalance,outT,outY] = Tenofovir_multipledose(p,OutputVar,TimeLen); 
%                 plot(outT,outY);
%                 title('Dosing Every 24 Hours for 10 Days: Case 3');
%                 xlabel('Time (hr)');
%                 ylabel('TFV in Central Compartment (nmol/L)');
                cmin=min(outY(end-100:end));
                cmax=max(outY(end-100:end));
                cminall(sc,dum)=cmin;
                cmaxall(sc,dum)=cmax;
            case 4
                ka= exp(-info(1,npeople)/120);
                Cl=Fbio*(number+136*info(2,dum)/info(3,dum));
                kout=info(4,dum);
                p = [dose V1 V2 Vcell Fbio ka Q Cl kout kcat1 kcat2 Km1 Km2 Km Vmax E1 E2 GSHi GSHe fu F Vbl N kleak]';
                [outAUC,outBalance,outT,outY] = Tenofovir_multipledose(p,OutputVar,TimeLen); 
%                 plot(outT,outY);
%                 title('Dosing Every 24 Hours for 10 Days: Case 4');
%                 xlabel('Time (hr)');
%                 ylabel('TFV in Central Compartment (nmol/L)');
                cmin=min(outY(end-100:end));
                cmax=max(outY(end-100:end));
                cminall(sc,dum)=cmin;
                cmaxall(sc,dum)=cmax;
        end
    end
end
%%
figure;
ax1=subplot(2,2,1);
hist(cmaxall(1,:));
title(ax1,'Cmax: Case 1')
xlabel(ax1,'TFV (nmol/L)')

ax1=subplot(2,2,2);
hist(cmaxall(2,:));
title(ax1,'Cmax: Case 2')
xlabel(ax1,'TFV (nmol/L)')

ax1=subplot(2,2,3);
hist(cmaxall(3,:));
title(ax1,'Cmax: Case 3')
xlabel(ax1,'TFV (nmol/L)')

ax1=subplot(2,2,4);
hist(cmaxall(4,:));
title(ax1,'Cmax: Case 4')
xlabel(ax1,'TFV (nmol/L)')