%% Driver for Part 4: Missed Dose Analysis
TimeLen = 15*24;
OutputVar = 1:15;
%% Missed Dose
dose = 300*10^6; %dose in ng, taken orally
p(1) = dose;
% metric - 1: AUC of TFV-DP, 2: Ctrough of TFV-DP, 3: Cmax of TFV-DP, 4: Steady State Viral Load
metric_missDose = zeros(4,TimeLen/24);
for missDose = 10
    [outMetric,Balance,t,y] = Tenofovir_missDose(p,p_viral,y0_viral,OutputVar,TimeLen,missDose);
    y(:,end+1) = 2*(y(:,14) + y(:,15))/(VD_virus*1000);
    y_missDose{missDose} = y;
%     tset{missDose} = t;
%     balanceset{missDose} = Balance;
    figure();
    plot(t,y(:,5)/(Vcell2*10^3),'linewidth',3)
    metric_missDose(:,missDose) = outMetric;
end

%% Retake Dose
dose = 300*10^6; %dose in ng, taken orally
p(1) = dose;
figure();
hold on;
metric_retakeDose = zeros(4,9);
for i = 1:9
    [outMetric,Balance,t,y] = Tenofovir_retakeDose(p,p_viral,y0_viral,OutputVar,TimeLen,10,i);
    y(:,end+1) = 2*(y(:,14) + y(:,15))/(VD_virus*1000);
    y_retakeDose{i} = y;
    plot(t,y(:,5)/(Vcell2*10^3),'linewidth',3)
    metric_retakeDose(:,i) = outMetric;
end
title('Concentration of TFV-DP in PBMC') %(zero = balance)
ylabel('TFV-DP (nmol/mL)')
xlabel('time (hrs)')
