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

%% Sensitivity Analysis
TimeLen=1/3;
p(1) = dose_set(2);
[metric0,~,~,~] = Tenofovir(p,p_viral,y0_viral,OutputVar,TimeLen);
metric0=metric0';
ParamDelta = 0.05;
metricsenspk = zeros(4,length(p));

%PK params
for i=1:length(p)
    psen=p;
    psen(i)=p(i)*(1.0+ParamDelta);
    [output,~,~,~] = Tenofovir(psen,p_viral,y0_viral,OutputVar,TimeLen);
    metricsenspk(:,i) = output; 
    Senspk(:,i)=((metricsenspk(:,i)-metric0)./metric0)/((psen(i)-p(i))/p(i)); % dim: metric by param
end
specs={'dose','V1','V2','Vcell','Fbio','ka','Q','Cl','kout','kcat1',...
    'kcat2','Km1','Km2','Km','Vmax','E1','E2','GSHi','GSHe','fu','Vbl','N','kleak'};
tt={'PK Parameters Local Sensitivity: AUC','PK Parameters Local Sensitivity: Ctrough',...
    'PK Parameters Local Sensitivity: Cmax','PK Parameters Local Sensitivity: Viral Load'};

for i=1:4
figure;
barh(Senspk(i,:));
xlim([-1,1]);
ylim([1,23]);
title(tt{i});
set(gca,'yticklabel',specs,'YTick', 1:1:23);
end

%PD params
metricsenspd = zeros(4,length(p_viral));
for i=1:length(p_viral)
    psen=p_viral;
    psen(i)=p_viral(i)*(1.0+ParamDelta);
    [output,~,~,~] = Tenofovir(p,psen,y0_viral,OutputVar,TimeLen);
    metricsenspd(:,i) = output; 
    Senspd(:,i)=((metricsenspd(:,i)-metric0)./metric0)/((psen(i)-p_viral(i))/p_viral(i)); % dim: metric by param
end

specs={'gammaT','gammaM','deltaT','deltaM','deltaPICT','deltaPICM','kT','kM','NhatT','NhatM',...
    'NT','NM','deltaT1','deltaT2','deltaM1','deltaM2','CL_n','CL_in','IC50'};
tt={'PD Parameters Local Sensitivity: AUC','PD Parameters Local Sensitivity: Ctrough',...
    'PD Parameters Local Sensitivity: Cmax','PD Parameters Local Sensitivity: Viral Load'};

for i=1:4
figure;
barh(Senspd(i,:));
xlim([-1,1]);
ylim([1,19]);
title(tt{i});
set(gca,'yticklabel',specs,'YTick', 1:1:19);
end
