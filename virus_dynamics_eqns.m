function dydt = virus_dynamics_eqns(t,y,p)

gammaT=p(1);
gammaM=p(2);
deltaT=p(3);
deltaM=p(4);
deltaPICT=p(5);
deltaPICM=p(6);
kT=p(7);
kM=p(8);
NhatT=p(9);
NhatM=p(10);
NT=p(11);
NM=p(12);
deltaT1=p(13);
deltaT2=p(14);
deltaM1=p(15);
deltaM2=p(16);
CL_n=p(17);
CL_in=p(18);


dydt=zeros(8,1);    % make it a column vector

% dydt(1-8), y(1-8) is in this order:
% TU: T cells
% MU: macrophages
% T1: infected T-cells prior to proviral genomic integration
% M1: infected macrophages prior to proviral genomic integration
% T2: infected T-cells after proviral genomic integration
% M2: infected macrophages after proviral genomic integration
% VI: free infectious virus
% VNI: free non-infectious virus

% beta=func_beta(etaterm);
% CL=func_CL(etaterm);
betaT=8*10^-12;
betaM=10^-14;
prev = 0.5;
CLT=((1/prev)-1)*betaT;
CLM=((1/prev)-1)*betaM;
 
 dydt(1) = gammaT+deltaPICT*y(3)-deltaT*y(1)-betaT*y(7)*y(1);
 dydt(2) = gammaM+deltaPICM*y(4)-deltaM*y(2)-betaM*y(7)*y(2);
 dydt(3) = betaT*y(7)*y(1)-(deltaT1+deltaPICT+kT)*y(3);
 dydt(4) = betaM*y(7)*y(2)-(deltaM1+deltaPICM+kM)*y(4);
 dydt(5) = kT*y(3)-deltaT2*y(5);
 dydt(6) = kM*y(4)-deltaM2*y(6);
 dydt(7) = NM*y(6)+NT*y(5)-y(7)*(CL_in+(CLT+betaT)*y(1)+(CLM+betaM)*y(2));
 dydt(8) = ((NhatT-NT)*y(5)+(NhatM-NM)*y(6))-CL_in*y(8);
end