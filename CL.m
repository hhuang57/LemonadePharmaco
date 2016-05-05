function CL = CL(etaterm)
% CLT, CLM: clearance of virus through unsuccessful infection of T-cells 
% and macrophages in the presence of TFV-DP
CL(1)=(1/0.5-etaterm)*8e-12; %T
CL(2)=(1/0.5-etaterm)*1e-14; %M