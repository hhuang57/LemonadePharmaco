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