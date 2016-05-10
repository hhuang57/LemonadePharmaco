%% Driver for Part 2
dose_set = [75, 150, 300, 600];
TimeLen = 15*24;
OutputVar = 1:15;

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

figure;
ax1=subplot(2,3,1);
semilogy(ax1,t,y(:,12),'linewidth',3)
title(ax1,'T2')
ylabel(ax1,'T2')
xlabel(ax1,'time (hrs)')

ax1=subplot(2,3,2);
semilogy(ax1,t,y(:,10),'linewidth',3)
title(ax1,'T1')
ylabel(ax1,'T1')
xlabel(ax1,'time (hrs)')

ax1=subplot(2,3,3);
semilogy(ax1,t,y(:,8),'linewidth',3)
title(ax1,'Tu')
ylabel(ax1,'Tu')
xlabel(ax1,'time (hrs)')

ax4=subplot(2,3,4);
semilogy(ax4,t,y(:,13),'linewidth',3)
title(ax4,'M2')
ylabel(ax4,'M2')
xlabel(ax4,'time (hrs)')

ax4=subplot(2,3,5);
semilogy(ax4,t,y(:,11),'linewidth',3)
title(ax4,'M1')
ylabel(ax4,'M1')
xlabel(ax4,'time (hrs)')

VD_virus = 50*3.1 + 9.6; % Volume of Distribution
ax4=subplot(2,3,6);
semilogy(ax4,t,2*(y(:,14)+y(:,15))/(VD_virus*1000),'linewidth',3)
title(ax4,'Viral Load')
ylabel(ax4,'# of HIV-1 RNA copies/mL')
xlabel(ax4,'time (hrs)')


% VD_virus = 50*3.1 + 9.6;
% ax4=subplot(2,3,6);
% plot(ax4,t,log10((y(:,14)+y(:,15)/(VD_virus*1000))/((y0_viral(7)+y0_viral(8))/(VD_virus*1000))),'linewidth',3)
% title(ax4,'Viral Load Decay')
% ylabel(ax4,'Log10 Viral Load Decay')
% xlabel(ax4,'time (hrs)')

% %% Step 4: Missed Dose Analysis
% %% Missed Dose
% dose = 300*10^6; %dose in ng, taken orally
% p(1) = dose;
% [outMetric,Balance,t,y] = Tenofovir_missDose(p,p_viral,y0_viral,OutputVar,TimeLen,3);
% AUC = outMetric(1)
% Ctrough = outMetric(2)
% Cmax = outMetric(3)
% figure();
% plot(t,y(:,5)/(Vcell2*10^3),'k','linewidth',3)
% title('Concentration of TFV-DP in PBMC') %(zero = balance)
% ylabel('TFV-DP(nmol/mL)')
% xlabel('time (hrs)')
% 
% %% Retake Dose
% dose = 300*10^6; %dose in ng, taken orally
% p(1) = dose;
% TimeLen = 8*24;
% figure();
% hold on;
% for i = 1:3
%     [outMetric,Balance,t,y] = Tenofovir_retakeDose(p,p_viral,y0_viral,OutputVar,TimeLen,3,i);
%     plot(t,y(:,5)/(Vcell2*10^3),'linewidth',3)
%     AUC(i) = outMetric(1);
%     Ctrough(i) = outMetric(2);
%     Cmax(i) = outMetric(3);
% end
% title('Concentration of TFV-DP in PBMC') %(zero = balance)
% ylabel('TFV-DP (nmol/mL)')
% xlabel('time (hrs)')
