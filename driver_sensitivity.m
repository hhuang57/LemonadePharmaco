%% driver for Sensitivity Analysis
TimeLen=15*24;
p(1) = 300*10^6;
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
