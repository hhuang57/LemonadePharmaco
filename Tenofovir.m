function [outMetric,outBalance,outT,outY] = Tenofovir(p,p_viral,y0_viral,OutputVar,TimeLen)

mtfv = 287.2; %Molecular weight of tenofovir (g/mol)
mmp = 367.2; %Molecular weight of tenofovir monophosphate
mdp = 447.18; %Molecular weight of tenofovir diphosphate
D0 = p(1)/mtfv; %nmol (oral)
V1 = p(2);
V2 = p(3);
Vcell = p(4);
Vbl = p(21);
N = p(22); %PBMC per L of blood
Vcell2 = Vcell*N*Vbl;
VD_virus = p_viral(20);
y0 = [0 0 0 0 0 D0 0]'; % moles
y0 = [y0; y0_viral'];
p = [p; p_viral];
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode23s(@Tenofovir_eqns,[0 24],y0,options,p);
DrugIn = ones(size(T1))*D0; % cumulative drug into system
TotalFreeD(:,1) = Y1(:,1);
TotalFreeD(:,2) = Y1(:,2);
TotalFreeD(:,3) = Y1(:,3);
TotalDmp = Y1(:,4);
TotalDdp = Y1(:,5);
TotalProD = Y1(:,6);
DrugOut = Y1(:,7); % cumulative drug eliminated from system
BalanceD1 = DrugIn - DrugOut - TotalProD - TotalFreeD(:,1) - TotalFreeD(:,2) - TotalFreeD(:,3) - TotalDmp - TotalDdp; %(zero = balance)
BalanceD1 = BalanceD1/max(DrugIn);
clearvars TotalFreeD;
for T=24:24:TimeLen
    y0 = Y1(end,:);
    y0(6) = y0(6) + D0;
    [t,y] = ode23s(@Tenofovir_eqns,[0 24],y0,options,p);
    DrugIn = ones(size(t))*(T/24 + 1)*D0; % cumulative drug into system
    TotalFreeD(:,1) = y(:,1);
    TotalFreeD(:,2) = y(:,2);
    TotalFreeD(:,3) = y(:,3);
    TotalDmp = y(:,4);
    TotalDdp = y(:,5);
    TotalProD = y(:,6);
    DrugOut = y(:,7); % cumulative drug eliminated from system
    balance = DrugIn - DrugOut - TotalProD - TotalFreeD(:,1) - TotalFreeD(:,2) - TotalFreeD(:,3) - TotalDmp - TotalDdp; %(zero = balance)
    balance = balance/max(DrugIn);
    BalanceD1 = [BalanceD1; balance];
    T1 = [T1; t+T];
    Y1 = [Y1; y];
    clearvars TotalFreeD;
end

% Include a check on the molecular balance. Don't want to do it visually
% for 48,000 runs! Instead define a criterion. For example, alert us for
% any mismatch by more than 1 molecule in a million (10^-6)
check = max(max(BalanceD1));
% fprintf ('Molecular Balance = %2.1e\n',check);
if check > 1.e-6
    fprintf ('*** Molecular Balance Violated ***\n');
end
outAUC=trapz(T1,Y1(:,5)/(Vcell2*10^3)); % AUC of TFV-DP
% DataInv = 1.01*max(Y1(:,5)) - Y1(:,5);
% [~,MinIdx] = findpeaks(DataInv);
% Ctrough = min(Y1(MinIdx,5))/(Vcell2*10^3);
Ctrough = min(Y1(T1 > (TimeLen - 24*5),5)/(Vcell2*10^3));
Cmax = max(Y1(:,5)/(Vcell2*10^3));
virLoad = 2*(Y1(end,14)+Y1(end,15))/(VD_virus*1000);
outMetric = [outAUC Ctrough Cmax virLoad];
outT=T1;
outY=Y1(:,OutputVar);
outBalance = BalanceD1;
end


