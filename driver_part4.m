%% Driver for Part 4: Missed Dose Analysis
load popCLIC50.mat
TimeLen = 15*24;
OutputVar = 1:15;
%% Missed Dose
dose = 300*10^6; %dose in ng, taken orally
p(1) = dose;
% metric - 1: AUC of TFV-DP, 2: Ctrough of TFV-DP, 3: Cmax of TFV-DP, 4: Steady State Viral Load
metric_missDose = zeros(4,length(popparam));
for j = 1:length(popparam)
    for missDose = 11
        p(8) = popparam(1,j);
        p_viral(19) = popparam(2,j);
        [outMetric,Balance,t,y] = Tenofovir_missDose(p,p_viral,y0_viral,OutputVar,TimeLen,missDose);
        metric_missDose(:,j) = outMetric;
    end
end
figure();
ax1=subplot(2,2,1);
hist(metric_missDose(1,:));
title(ax1,'AUC of TFV-DP for 11th Missed Dose')
ylabel(ax1,'Incidence')
xlabel(ax1,'AUC of TFV-DP (nmol*hr/mL)');

ax2=subplot(2,2,2);
hist(metric_missDose(2,:))
title(ax2,'Ctrough of TFV-DP for 11th Missed Dose') %(zero = balance)
ylabel(ax2,'Incidence')
xlabel(ax2,'Ctrough of TFV-DP (nmol/mL)')

ax3 = subplot(2,2,3);
hist(metric_missDose(3,:))
title(ax3,'Cmax of TFV-DP for 11th Missed dose')
ylabel(ax3,'Incidence')
xlabel(ax3,'Cmax of TF-DP (nmol/L)')

ax4=subplot(2,2,4);
hist(metric_missDose(4,:))
title(ax4,'Steady State Viral Load for 11th Missed Dose') %(zero = balance)
ylabel(ax4,'Incidence')
xlabel(ax4,'# of HIV-1 RNA copies/mL')
%% Retake Dose
dose = 300*10^6; %dose in ng, taken orally
p(1) = dose;
figure();
hold on;
metric_retakeDose = zeros(4,9);
for i = 1:9
    [outMetric,Balance,t,y] = Tenofovir_retakeDose(p,p_viral,y0_viral,OutputVar,TimeLen,11,i);
    y(:,end+1) = 2*(y(:,14) + y(:,15))/(VD_virus*1000);
    y_retakeDose{i} = y;
    plot(t/24,y(:,5)/(Vcell2*10^3),'linewidth',3)
    metric_retakeDose(:,i) = outMetric;
end
title('11th Missed Dose and Retaken Dose Results') %(zero = balance)
ylabel('TFV-DP in PBMC (nmol/mL)')
xlabel('time (days)')

save driverpart4.mat metric_missDose metric_retakeDose;