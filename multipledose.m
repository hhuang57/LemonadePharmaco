function [outAUC,outBalance,outT,outY] = Tenofovirpop(p,OutputVar,TimeLen)

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

y0 = [0 0 0 0 0 D0 0]'; % moles
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
% [T1,Y1] = ode45(@Tenofovir_eqns,[0 TimeLen],y0,options,p);
[T1,Y1] = ode23s(@Tenofovir_eqns_pop,[0 24],y0,options,p);
% DrugIn = ones(size(T1))*D0; % cumulative drug into system
% TotalFreeD(:,1) = Y1(:,1);
% TotalFreeD(:,2) = Y1(:,2);
% TotalFreeD(:,3) = Y1(:,3);
% TotalDmp = Y1(:,4);
% TotalDdp = Y1(:,5);
% TotalProD = Y1(:,6);
% DrugOut = Y1(:,7); % cumulative drug eliminated from system
% BalanceD1 = DrugIn - DrugOut - TotalProD - TotalFreeD(:,1) - TotalFreeD(:,2) - TotalFreeD(:,3) - TotalDmp - TotalDdp; %(zero = balance)
% % Include a check on the molecular balance. Don't want to do it visually
% % for 48,000 runs! Instead define a criterion. For example, alert us for
% % any mismatch by more than 1 molecule in a million (10^-6)
% check = max(max(BalanceD1))/(D0/V1);
% % fprintf ('Molecular Balance = %2.1e\n',check);
% if check > 1.e-5
%     fprintf ('*** Molecular Balance Violated ***\n');
%end

y0_2 = [Y1(end,1) Y1(end,2) Y1(end,3) Y1(end,4) Y1(end,5) Y1(end,6)+D0 Y1(end,7)]'; % moles
[T1_2,Y1_2] = ode23s(@Tenofovir_eqns_pop,[24 48],y0_2,options,p);
Y1 = [Y1; Y1_2];
T1 = [T1; T1_2];
% TotalD1(:,1) = Y1_2(:,1)*Vp;
% TotalD1(:,2) = Y1_2(:,2)*Vm;
% TotalD1(:,3) = Y1_2(:,3)*Vb;
% TotalD1(:,4) = Y1_2(:,4)*Vp;
% freeD1 = TotalD1(:,1)+TotalD1(:,2)+TotalD1(:,3);
% DrugIn1 =  0; % cumulative drug into system
% DrugOut1 = Y1_2(:,6)*Vcl ; % cumulative drug eliminated from system
% BalanceD1_2 = DrugIn1 - DrugOut1 - TotalD1(:,1) - TotalD1(:,2) - TotalD1(:,3)- TotalD1(:,4)+ y0_2(1)*Vp + y0_2(2)*Vm + y0_2(3)*Vb + y0_2(4)*Vp + y0_2(6)*Vcl; %(zero = balance)
% BalanceD1 = [BalanceD1; BalanceD1_2];

y0_3 = [Y1_2(end,1) Y1_2(end,2) Y1_2(end,3) Y1_2(end,4) Y1_2(end,5) Y1_2(end,6)+D0 Y1_2(end,7)]'; % moles
[T1_3,Y1_3] = ode23s(@Tenofovir_eqns_pop,[48 72],y0_3,options,p);
Y1 = [Y1; Y1_3];
T1 = [T1; T1_3];
% TotalD1(:,1) = Y1_3(:,1)*Vp;
% TotalD1(:,2) = Y1_3(:,2)*Vm;
% TotalD1(:,3) = Y1_3(:,3)*Vb;
% TotalD1(:,4) = Y1_3(:,4)*Vp;
% freeD1 = TotalD1(:,1)+TotalD1(:,2)+TotalD1(:,3);
% DrugIn1 =  0; % cumulative drug into system
% DrugOut1 = Y1_3(:,6)*Vcl ; % cumulative drug eliminated from system
% BalanceD1_3 = DrugIn1 - DrugOut1 - TotalD1(:,1) - TotalD1(:,2) - TotalD1(:,3)- TotalD1(:,4)+ y0_2(1)*Vp + y0_2(2)*Vm + y0_2(3)*Vb + y0_2(4)*Vp + y0_2(6)*Vcl; %(zero = balance)
% BalanceD1 = [BalanceD1; BalanceD1_3];

outAUC=trapz(T1,Y1(:,1));
outT=T1;
outY=Y1(:,1);
outBalance = BalanceD1;
end
