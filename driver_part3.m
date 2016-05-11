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
LPVdistrib=[0.31,0.69];

npeople=10;
info=zeros(5,npeople); %age,weight,scr,kout,Cl

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

for p5=1:npeople
    a=rand();
    if a<LPVdistrib(1)
        info(5,p5)=0.75*p(8);
    else
        info(5,p5)=p(8);
    end
end

figure;
axx1=subplot(2,3,1);
hist(info(1,:));
title('Age Distribution');
axx2=subplot(2,3,2);
hist(info(2,:));
title('Weight Distribution');
axx3=subplot(2,3,3);
hist(info(3,:));
title('SCR Distribution');
axx4=subplot(2,3,4);
hist(info(4,:));
title('kout Distribution');
axx5=subplot(2,3,5);
hist(info(5,:));
title('Cl Distribution');

%% PopPK
dose_set = [75, 150, 300, 600]*10^6; % Dose in mg converted to ng
TimeLen = 15*24;
OutputVar = 1:15;
p(1) = dose_set(2);
[metric0,~,~,~] = Tenofovir(p,p_viral,y0_viral,OutputVar,TimeLen);

auc = zeros(5,npeople); %base line, vary ka, vary cl, vary kout, vary cl
auc(1,:)=repmat(metric0(1),[1,npeople]);
ctro = zeros(5,npeople);
ctro(1,:)=repmat(metric0(2),[1,npeople]);
cmax = zeros(5,npeople);
cmax(1,:)=repmat(metric0(3),[1,npeople]);
vload = zeros(5,npeople);
vload(1,:)=repmat(metric0(4),[1,npeople]);

for sc=1:4
    for dum=1:npeople
        switch sc
            case 1
                p1=p;
                p1(6)= exp(-info(1,dum)/50);
                [output,Balance,t,y] = Tenofovir(p1,p_viral,y0_viral,OutputVar,TimeLen); 
                auc(2,dum)=output(1);
                ctro(2,dum)=output(2);
                cmax(2,dum)=output(3);
                vload(2,dum)=output(4);
                
            case 2
                p2=p;
                p2(8)=Fbio*(number+136*info(2,dum)/info(3,dum));
                [output,Balance,t,y] = Tenofovir(p2,p_viral,y0_viral,OutputVar,TimeLen); 
                auc(3,dum)=output(1);
                ctro(3,dum)=output(2);
                cmax(3,dum)=output(3);
                vload(3,dum)=output(4);
                
            case 3
                p3=p;
                p3(9)=info(4,dum);
                [output,Balance,t,y] = Tenofovir(p3,p_viral,y0_viral,OutputVar,TimeLen); 
                auc(4,dum)=output(1);
                ctro(4,dum)=output(2);
                cmax(4,dum)=output(3);
                vload(4,dum)=output(4);
                
            case 4
                p4=p;
                p4(8)=info(5,dum);
                [output,Balance,t,y] = Tenofovir(p4,p_viral,y0_viral,OutputVar,TimeLen); 
                auc(5,dum)=output(1);
                ctro(5,dum)=output(2);
                cmax(5,dum)=output(3);
                vload(5,dum)=output(4);
        end
    end
end

%%
specs={'Baseline','ka','CL_1','kout','CL_2'};
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
