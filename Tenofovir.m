function out = Tenofovir(?)

TimeLen = 24; % hours
mtfv = 287.2*10^3; %Molecular weight of tenofovir (mg/mol)
mmp = 367.2*10^3; %Molecular weight of tenofovir monophosphate
mdp = 447.18*10^3; %Molecular weight of tenofovir diphosphate

D0 = 300; %mg (oral)
y0 = [D0/VD 0]'; % mg/kg*L
    % 1 = drug in blood; 2 = drug in degr
p = [CL VD]'; % parameter array

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode45(@Artesumate_eqns,[0 TimeLen],y0,options,p);

DrugInBolus = ones(length(T1),1)*(D0/VD) ;
TotalFreeD(:,1) = Y1(:,1)*VD;
DrugIn = DrugInBolus*VD ; % cumulative drug into system
DrugOut = Y1(:,2)*VD ; % cumulative drug eliminated from system
BalanceD1 = DrugIn - DrugOut - TotalFreeD(:,1); %(zero = balance)

% Include a check on the molecular balance. Don't want to do it visually
% for 48,000 runs! Instead define a criterion. For example, alert us for 
% any mismatch by more than 1 molecule in a million (10^-6)
check = max(max(BalanceD1))/(D0/VD);
% fprintf ('Molecular Balance = %2.1e\n',check);
if check > 1.e-6
    fprintf ('*** Molecular Balance Violated ***\n');
end

end



