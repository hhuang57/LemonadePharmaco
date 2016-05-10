%% Driver for Part 2
dose_set = [75, 150, 300, 600]*10^6; % Dose in mg converted to ng
TimeLen = 15*24;
OutputVar = 1:15;
% metric - 1: AUC of TFV-DP, 2: Ctrough of TFV-DP, 3: Cmax of TFV-DP, 4: Steady State Viral Load
metric = zeros(4,length(dose_set));
for i = 1:4
    dose = dose_set(i); % taken orally
    p(1) = dose;
    [output,Balance,t,y] = Tenofovir(p,p_viral,y0_viral,OutputVar,TimeLen);
    y(:,end+1) = 2*(y(:,14) + y(:,15))/(VD_virus*1000);
    yset{i} = y;
    tset{i} = t/24; %Convert time from hours to days
    metric(:,i) = output; 
    balanceset{i} = Balance;
end
y1 = yset{1};
y2 = yset{2};
y3 = yset{3};
y4 = yset{4};

%% Plot Time-Dependent Variables
figure;
ax1=subplot(2,2,1);
plot(ax1,tset{1},y1(:,1)/(V1*10^3),'k',tset{2},y2(:,1)/(V1*10^3),'r',tset{3},y3(:,1)/(V1*10^3),'b',tset{4},y4(:,1)/(V1*10^3),'g','linewidth',3)
title(ax1,'Concentration of TFV in Central Compartment')
ylabel(ax1,'TFV (nmol/mL)')
xlabel(ax1,'time (days)');

ax2=subplot(2,2,2);
plot(ax2,tset{1},y1(:,5)/(Vcell2*10^3),'k',tset{2},y2(:,5)/(Vcell2*10^3),'r',tset{3},y3(:,5)/(Vcell2*10^3),'b',tset{4},y4(:,5)/(Vcell2*10^3),'g','linewidth',3)
title(ax2,'Concentration of TFV-DP in PBMC') %(zero = balance)
ylabel(ax2,'TFV-DP (nmol/mL)')
xlabel(ax2,'time (days)')

ax3 = subplot(2,2,3);
semilogy(ax3,tset{1},y1(:,16),'k',tset{2},y2(:,16),'r',tset{3},y3(:,16),'b',tset{4},y4(:,16),'g','linewidth',3)
title(ax3,'Viral Load Decay')
ylabel(ax3,'# of HIV-1 RNA copies/mL')
xlabel(ax3,'time (days')

ax4=subplot(2,2,4);
plot(ax4,tset{1},balanceset{1},'k',tset{2},balanceset{2},'r',tset{3},balanceset{3},'b',tset{4},balanceset{4},'g','linewidth',3)
title(ax4,'Molecular Balance') %(zero = balance)
ylabel(ax4,'Balance of Drug (nmol)')
xlabel(ax4,'time (days)')

legend('75 mg','150 mg','300 mg','600 mg');

save driver2.mat tset yset metric;
