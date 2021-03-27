%% Setup Parmeters
clear all; close all;

f = 27.12e6;    %frequency of operation
w = 2*pi*f;
Rl = 500;        % Load on the secondary coil (<37.5)
Ql = 10;        % Quality Factor for class E amplifier
K = 0.2;          % Coupling Coefficient
Ciss = 14.5e-12;% EPC2037: 14.5e-12; EPC8002: 20.5e-12
Coss = 17e-12;  % EPC2037: 17e-12; EPC8002: 13e-12
Crss = 1e-12;   % EPC2037: 1e-12; EPC8002: 1e-12
Co = Coss - Crss;% Typical Output Capacitance of the Switch

Lp = 380e-9;   % Primary coil Inductance
Rp = 0.32;      % Primary coil parasitic series Resistance
Cpp = 1.45e-12; % Primary coil parasitic parallel Capacitance

Ls = 240e-9;   % Secondary coil Inductance
Rs = 0.28;      % Secondary coil parasitic series Resistance
Cps = 1.26e-12; % Secondary coil parasitic parallel Capacitance
Cs = 144e-12;   % Capacitor for parallel resonance on the secondary side

Lc = 10e-6;     % Choke Inductance
Lcsi = 0.2e-9;  % Common Source Inductance
Rg = 1;         % Switch Gate Resistance

%% Transformer T - model

M = K*(Ls*Lp)^0.5;  % Mutual Inductance
L1 = Lp - M;
L2 = Ls - M;
Z1 = Rs + j*w*L2 + Rl/(1 + j*w*Rl*(Cps + Cs));
Z2 = Rp + j*w*L1 + j*w*M*Z1/(Z1 + j*w*M);
Z = Z2/(1 + j*w*Z2*Cpp)

%% Calculate CMCD parameters
% C = Ql/(w*real(Z)) - Co/2 + 8e-12
% L = real(Z)/(w*Ql) - imag(Z)/w

C = Ql/(w*20) - Co/2 + 8e-12
L = 20/(w*Ql)

%% Writing into the file
fileID = fopen('CMCD_param.txt','w');

% fprintf(fileID,'.tran 0 36.8u 27.6u\n');
% fprintf(fileID,';ac dec 4 1K 1G\n');
fprintf(fileID,'.four %dHz I(R1)\n\n',f);

fprintf(fileID,'.param f=%d Rl=%d Ql=%d\n', f,Rl,Ql);
fprintf(fileID,'.param L=%d C=%d\n',L,C);        % Class E Param
fprintf(fileID,'.param Lp=%d Rp=%d Cpp=%d\n',Lp,Rp,Cpp);  % Primary Coil Param
fprintf(fileID,'.param Ls=%d Rs=%d Cps=%d\n',Ls,Rs,Cps);  % Secondary Coil Param
fprintf(fileID,'.param Lc=%d Lcsi=%d Rg=%d Cs=%d\n',Lc,Lcsi,Rg,Cs);% Other Param
%fprintf(fileID,'K Lp Ls %d\n',K);
fprintf(fileID,';.step param C 1.1e-10 1.2e-10 0.02e-10\n\n');

fprintf(fileID,'.meas Pout AVG (I(R1)*(V(Vo+)-V(Vo-)))\n');
fprintf(fileID,'.meas Pin AVG (I(V1)*V(Vdd))\n');
fprintf(fileID,'.meas Pgate AVG (I(V2)*V(Vg+) + I(V3)*V(Vg-))\n');
fprintf(fileID,'.meas eff param Pout/Pin\n');
fprintf(fileID,'.meas PAE param Pout/(Pin+Pgate)\n');
fclose(fileID);
