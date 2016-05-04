function dydt = Tenofovir_eqns(t,y,p)
mtfv = 287.2*10^3; %Molecular weight of tenofovir (mg/mol)
mmp = 367.2*10^3; %Molecular weight of tenofovir monophosphate
mdp = 447.18*10^3; %Molecular weight of tenofovir diphosphate
D0=p(1);
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
N = p(23); %PBMC per ml of blood

kcl = Cl/V1; %clearance rate of TFV
k12 = Q/V1; %transfer rate constants from central to peripheral
k21 = Q/V2; %transfer rate constants from peripheral to central

dydt=zeros(7,1);    % make it a column vector (e.g. (3,1)

% Units are in concentration (mg/L)
% 1 = TFV in central compartment
% 2 = TFV in peripheral compartment
% 3 = TFV in PBMC
% 4 = TFV-MP in PBMC
% 5 = TFV-DP in PBMC
% 6 = TDF in dosing compartment
% 7 = clearance???

 dTFVin = [(Vmax*fu*(1-F)*y(1))/(Km + fu*(1-F)*y(1))]*((GSHi-GSHe)/GSHe)*(V1/(Vcell*N*Vbl));
 dTFVef = (Vmax*fu*(1-F)*y(3))/(Km + fu*(1-F)*y(3)/(mtfv*Vcell*N*Vbl))*(Vcell*N*Vbl/V1);
 kmp = (kcat1*E1)/(Km1 + y(3));
 kdp = (kcat2*E2)/(Km2 + y(4));
 
 dydt(1) = Fbio*ka*y(6) - y(1)*kcl - y(1)*k12 + y(2)*k21*(V2/V1) +...
     y(3)*kleak*(Vcell*N*Vbl/V1) - dTFVin + dTFVef;
 dydt(2) = y(1)*k12*(V1/V2) - y(2)*k21;
 dydt(3) = - y(3)*kleak + dTFVin - dTFVef - y(3)*kmp;
 dydt(4) = y(3)*kmp - y(4)*kdp;
 dydt(5) = y(4)*kdp - y(4)*kout;
 dydt(6) = -y(6)*ka;
 dydt(7) = y(1)*kcl + y(6)*(1-Fbio)*ka + y(4)*kout*(Vcell*N*Vbl/V1);
 
%  dTFVin = [(Vmax*fu*(1-F)*y(1))/(Km + fu*(1-F)*y(1)/(mtfv*V1))]*((GSHi-GSHe)/GSHe);
%  dTFVef = (Vmax*fu*(1-F)*y(3))/(Km + fu*(1-F)*y(3)/(mtfv*Vcell*N*Vbl));
%  kmp = (kcat1*E1)/(Km1 + y(3)/(mtfv*Vcell*N*Vbl));
%  kdp = (kcat2*E2)/(Km2 + y(4)/(mmp*Vcell*N*Vbl));
%  
%  dydt(1) = Fbio*ka*y(6) - y(1)*kcl - y(1)*k12 + y(2)*k21 +...
%      y(3)*kleak - dTFVin + dTFVef;
%  dydt(2) = y(1)*k12 - y(2)*k21;
%  dydt(3) = - y(3)*kleak + dTFVin - dTFVef - y(3)*kmp;
%  dydt(4) = y(3)*kmp - y(4)*kdp;
%  dydt(5) = y(4)*kdp - y(4)*kout;
%  dydt(6) = -y(6)*ka;
%  dydt(7) = y(1)*kcl + y(6)*(1-Fbio)*ka + y(4)*kout
end

