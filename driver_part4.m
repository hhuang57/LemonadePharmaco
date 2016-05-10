%% Driver for Part 4: Missed Dose Analysis

%% Missed Dose
dose = 300*10^6; %dose in ng, taken orally
p(1) = dose;
% metric - 1: AUC of TFV-DP, 2: Ctrough of TFV-DP, 3: Cmax of TFV-DP, 4: Steady State Viral Load
metric = zeros(4,TimeLen/24);
for missDose = 1:(TimeLen/24)
    [outMetric,Balance,t,y] = Tenofovir_missDose(p,p_viral,y0_viral,OutputVar,TimeLen,3);
    y(:,end+1) = 2*(y(:,14) + y(:,15))/(VD_virus*1000);
    yset{missDose} = y;
    tset{missDose} = t;
    balanceset{missDose} = Balance;
    metric(:,missDose) = outMetric;
end
% figure();
% plot(t,y(:,5)/(Vcell2*10^3),'k','linewidth',3)
% title('Concentration of TFV-DP in PBMC') %(zero = balance)
% ylabel('TFV-DP(nmol/mL)')
% xlabel('time (hrs)')

%% Retake Dose
dose = 300*10^6; %dose in ng, taken orally
p(1) = dose;
TimeLen = 8*24;
figure();
hold on;
for i = 1:3
    [outMetric,Balance,t,y] = Tenofovir_retakeDose(p,p_viral,y0_viral,OutputVar,TimeLen,3,i);
    plot(t,y(:,5)/(Vcell2*10^3),'linewidth',3)
    AUC(i) = outMetric(1);
    Ctrough(i) = outMetric(2);
    Cmax(i) = outMetric(3);
end
title('Concentration of TFV-DP in PBMC') %(zero = balance)
ylabel('TFV-DP (nmol/mL)')
xlabel('time (hrs)')
