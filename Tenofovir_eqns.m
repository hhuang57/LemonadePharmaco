function dydt = Tenofovir_eqns(t,y,p)

V1=p(2); %Volume of central compartment
V2=p(3); %Volume of peripheral compartment
Vcell=p(4); %Volume of PBMC compartment
Fbio = p(5); %Bioavailability of tenofovir
ka = p(6); %Absorption constant
Q = p(7); %Intercompartment clearance (L/h)
Cl = p(8); %Clearance
kout = p(9); %Elimination of TFV-DP from PMBC
kcat1 = p(10); %Turnover number for AK2 and TFV
kcat2 = p(11); %Turnover number for NDKA and TFV-MP
Km1 = p(12); %Michaelis-Menten constant for AK2 and TFV
Km2 = p(13); %Michaelis-Menten constant for NDKA and TFV-MP
Km = p(14); %Michaelis-Menten constant for TFV uptake
Vmax = p(15); %Maximum rate of TFV uptake
E1 = p(16); %Enzyme concentration for hAK2
E2 = p(17); %Enzyme concentration for NDA
GSHi = p(18); %Intracellular GSH concentration
GSHe = p(19); %Extracellular GSH concentration
fu = p(20); %unbound fraction of TFV
F = p(21); %unionized fraction of TFV
Vbl = p(22);
N = p(23); %PBMC per L of blood
kleak = p(24); % Leakage of TFV
kcl = Cl/V1; %clearance rate of TFV
k12 = Q/V1; %transfer rate constants from central to peripheral
k21 = Q/V2; %transfer rate constants from peripheral to central
Vcell2 = Vcell*N*Vbl;

dydt=zeros(7,1);    % make it a column vector (e.g. (3,1)

% Units are in concentration (nmol/L)
% 1 = TFV in central compartment
% 2 = TFV in peripheral compartment
% 3 = TFV in PBMC
% 4 = TFV-MP in PBMC
% 5 = TFV-DP in PBMC
% 6 = TDF in dosing compartment
% 7 = clearance???

 dTFVin = [(Vmax*fu*y(1))/(Km + fu*y(1)/V1)]*((GSHi-GSHe)/GSHe);
 dTFVef = (Vmax*fu*y(3))/(Km + fu*y(3)/Vcell2);
 kmp = (kcat1*E1)/(Km1 + (y(3)/Vcell2));
 kdp = (kcat2*E2)/(Km2 + (y(4)/Vcell2));
 
 dydt(1) = Fbio*ka*y(6) - y(1)*kcl - y(1)*k12 + y(2)*k21 +...
     y(3)*kleak - dTFVin + dTFVef;
 dydt(2) = y(1)*k12 - y(2)*k21;
 dydt(3) = - y(3)*kleak + dTFVin - dTFVef - y(3)*kmp;
 dydt(4) = y(3)*kmp - y(4)*kdp;
 dydt(5) = y(4)*kdp - y(5)*kout;
 dydt(6) = -y(6)*Fbio*ka;
 dydt(7) = y(1)*kcl + y(5)*kout;
end

