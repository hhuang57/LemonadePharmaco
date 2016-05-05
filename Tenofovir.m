function [outAUC,outBalance,outT,outY] = Tenofovir(p,OutputVar,TimeLen)

mtfv = 287.2; %Molecular weight of tenofovir (g/mol)
mmp = 367.2; %Molecular weight of tenofovir monophosphate
mdp = 447.18; %Molecular weight of tenofovir diphosphate
D0 = p(1)/mtfv; %nmol (oral)
V1 = p(2);
V2 = p(3);
Vcell = p(4);
Vbl = p(22);
N = p(23); %PBMC per L of blood
Vcell2 = Vcell*N*Vbl;

% y0 = [0 0 0 0 0 D0/V1 0]'; % moles/L
y0 = [0 0 0 0 0 D0 0]'; % moles 
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
% [T1,Y1] = ode45(@Tenofovir_eqns,[0 TimeLen],y0,options,p);
[T1,Y1] = ode23s(@Tenofovir_eqns,[0 TimeLen],y0,options,p);

DrugIn = ones(size(T1))*D0; % cumulative drug into system
TotalFreeD(:,1) = Y1(:,1);
TotalFreeD(:,2) = Y1(:,2);
TotalFreeD(:,3) = Y1(:,3);
TotalDmp = Y1(:,4);
TotalDdp = Y1(:,5);
TotalProD = Y1(:,6);
DrugOut = Y1(:,7); % cumulative drug eliminated from system
BalanceD1 = DrugIn - DrugOut - TotalProD - TotalFreeD(:,1) - TotalFreeD(:,2) - TotalFreeD(:,3) - TotalDmp - TotalDdp; %(zero = balance)

% Include a check on the molecular balance. Don't want to do it visually
% for 48,000 runs! Instead define a criterion. For example, alert us for 
% any mismatch by more than 1 molecule in a million (10^-6)
check = max(max(BalanceD1))/(D0/V1);
% fprintf ('Molecular Balance = %2.1e\n',check);
if check > 1.e-6
    fprintf ('*** Molecular Balance Violated ***\n');
end

outAUC=trapz(T1,Y1(:,OutputVar));
outT=T1;
outY=Y1(:,OutputVar);
outBalance = BalanceD1;
end


