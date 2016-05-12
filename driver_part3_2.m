%% Population

%http://www.cdc.gov/hiv/statistics/overview/
agedistrib=[0.1054,0.1957,0.2269,0.383,0.089]; %15-19,25-29,35-39,45-49,55-59
cumagedistrib=[0.089,0.383+0.089,0.2269+0.383+0.089,...
    0.1957+0.2269+0.383+0.089,0.1054+0.1957+0.2269+0.383+0.089]; 

%things for CL
BWstdev = [25.46,40.82,40.82,29.45,29.45]*0.454; %kg
BWmean = [96.8,169.6,169.6,173.2,173.2]*0.454; %kg

scrstdev=19; %muM
scrmean=83; %muM
number=36.2;

lb=.04*10^3;
ub=8.5*10^3;
icmean=(ub+lb)/2;
icstdev=(ub-lb)/6; % ub and lb both 3 stdev away from the mean


npeople=10; %100
info=zeros(4,npeople); %age,weight,scr,IC50

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
    info(4,p4)=icstdev.*randn(1,1)+icmean;
end

figure;
axx1=subplot(2,2,1);
hist(info(1,:));
title('Age Distribution');
specs={'Age 15-19','Age 25-29','Age 35-39','Age 45-49','Age 55-59'};
set(gca,'xticklabel',specs,'XTick', 17:10:57);
axx2=subplot(2,2,2);
hist(info(2,:));
title('Weight Distribution');
xlabel('Weight (kg)');
axx3=subplot(2,2,3);
hist(info(3,:));
title('SCR Distribution');
xlabel('SCR Concentration (muM)');
axx4=subplot(2,2,4);
hist(info(4,:));
title('IC50 Distribution');
xlabel('IC50 (nmol/L)');

popparam=zeros(2,npeople);
popparam(1,:)=Fbio*(number+136/3*info(2,:)./info(3,:)); %CL in PK, added /3 to literature eqn to account for discrepancy in lit values
popparam(2,:)=info(4,:); %IC50
save popCLIC50.mat popparam;
%%
dose_set = [75, 150, 300, 600]*10^6; % Dose in mg converted to ng
TimeLen = 15*24;
OutputVar = 1:15;
p(1) = dose_set(2);
p(8)= median(popparam(1,:)); 
p_viral(19)=icmean;
[metric0,~,~,~] = Tenofovir(p,p_viral,y0_viral,OutputVar,TimeLen);

auc = zeros(4,npeople); %base line, vary IC50, vary cl, vary both
auc(1,:)=repmat(metric0(1),[1,npeople]);
ctro = zeros(4,npeople);
ctro(1,:)=repmat(metric0(2),[1,npeople]);
cmax = zeros(4,npeople);
cmax(1,:)=repmat(metric0(3),[1,npeople]);
vload = zeros(4,npeople);
vload(1,:)=repmat(metric0(4),[1,npeople]);

for sc=1:4
    for dum=1:npeople
        switch sc
            case 1
                p1=p_viral;
                p1(19)= popparam(2,dum);
                [output,Balance,t,y] = Tenofovir(p,p1,y0_viral,OutputVar,TimeLen); 
                auc(2,dum)=output(1);
                ctro(2,dum)=output(2);
                cmax(2,dum)=output(3);
                vload(2,dum)=output(4);
                
            case 2
                p2=p;
                p2(8)=popparam(1,dum); 
                [output,Balance,t,y] = Tenofovir(p2,p_viral,y0_viral,OutputVar,TimeLen); 
                auc(3,dum)=output(1);
                ctro(3,dum)=output(2);
                cmax(3,dum)=output(3);
                vload(3,dum)=output(4);
                
            case 3
                p3=p_viral;
                p4=p;
                p3(19)=popparam(2,dum);
                p4(8)=popparam(1,dum);
                [output,Balance,t,y] = Tenofovir(p4,p3,y0_viral,OutputVar,TimeLen); 
                auc(4,dum)=output(1);
                ctro(4,dum)=output(2);
                cmax(4,dum)=output(3);
                vload(4,dum)=output(4);
                
        end
    end
end

%%
specs={'Baseline','IC50','CL','Both'};
figure;
boxplot(auc'); %one column is a group
title('AUC of TFV-DP');
ylabel('hr-nmol/mL');
set(gca,'xticklabel',specs,'XTick', 1:1:5);
figure;
boxplot(ctro');
title('Ctrough of TFV-DP');
ylabel('nmol/mL');
set(gca,'xticklabel',specs,'XTick', 1:1:5);
figure;
boxplot(cmax');
title('Cmax of TFV-DP');
ylabel('nmol/mL');
set(gca,'xticklabel',specs,'XTick', 1:1:5);
figure;
boxplot(vload');
title('Viral Load');
ylabel('# of HIV-1 RNA copies/mL');
set(gca,'xticklabel',specs,'XTick', 1:1:5);
