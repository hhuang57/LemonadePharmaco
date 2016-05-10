function dydt = Tenofovir_eqns(t,y,p)
p_viral = p(23+1:end);
y_viral = y(8:end);
% PK component
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
Vbl = p(21);
N = p(22); %PBMC per L of blood
kleak = p(23); % Leakage of TFV
kcl = Cl/V1; %clearance rate of TFV
k12 = Q/V1; %transfer rate constants from central to peripheral
k21 = Q/V2; %transfer rate constants from peripheral to central
Vcell2 = Vcell*N*Vbl;

dydt=zeros(7+8,1);    % make it a column vector (e.g. (3,1)

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
 
% PD component
gammaT=p_viral(1);
gammaM=p_viral(2);
deltaT=p_viral(3);
deltaM=p_viral(4);
deltaPICT=p_viral(5);
deltaPICM=p_viral(6);
kT=p_viral(7);
kM=p_viral(8);
NhatT=p_viral(9);
NhatM=p_viral(10);
NT=p_viral(11);
NM=p_viral(12);
deltaT1=p_viral(13);
deltaT2=p_viral(14);
deltaM1=p_viral(15);
deltaM2=p_viral(16);
CL_n=p_viral(17);
CL_in=p_viral(18);
IC50 = p_viral(19);
prev=0.5;
reducTerm = IC50/(IC50 + (y(5)/Vcell2)); % 1 - n(t)
betaT0 = 8*10^-12;
betaM0 = 10^-14;
betaT = betaT0*reducTerm;
betaM = betaM0*reducTerm;
CLT = ((1/prev)-reducTerm)*betaT;
CLM = ((1/prev)-reducTerm)*betaM;

% 1,8: TU: T cells
% 2,9: MU: macrophages
% 3,10: T1: infected T-cells prior to proviral genomic integration
% 4,11: M1: infected macrophages prior to proviral genomic integration
% 5,12: T2: infected T-cells after proviral genomic integration
% 6,13: M2: infected macrophages after proviral genomic integration
% 7,14: VI: free infectious virus
% 8,15: VNI: free non-infectious virus
 dydt(1+7) = gammaT+deltaPICT*y_viral(3)-deltaT*y_viral(1)-betaT*y_viral(7)*y_viral(1);
 dydt(2+7) = gammaM+deltaPICM*y_viral(4)-deltaM*y_viral(2)-betaM*y_viral(7)*y_viral(2);
 dydt(3+7) = betaT*y_viral(7)*y_viral(1)-(deltaT1+deltaPICT+kT)*y_viral(3);
 dydt(4+7) = betaM*y_viral(7)*y_viral(2)-(deltaM1+deltaPICM+kM)*y_viral(4);
 dydt(5+7) = kT*y_viral(3)-deltaT2*y_viral(5);
 dydt(6+7) = kM*y_viral(4)-deltaM2*y_viral(6);
 dydt(7+7) = NM*y_viral(6)+NT*y_viral(5)-y_viral(7)*(CL_in+(CLT+betaT)*y_viral(1)+(CLM+betaM)*y_viral(2));
 dydt(8+7) = ((NhatT-NT)*y_viral(5)+(NhatM-NM)*y_viral(6))-CL_in*y_viral(8);

end

