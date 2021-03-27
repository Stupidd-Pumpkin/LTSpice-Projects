%% Instructions
% All the parameters listed below can be tweaked for finding the required component values.
% Please be careful to check the impedence values Z1, Z2, C1, C etc are legible.
% Change the Quality factor to as per need if a bottleneck is reached in the capacitor and inductor values.
% The cable properties correspond to the cable used by the HCL team.
% Check the output file before performing the LTSpice simulations based on the MATLAB results.
%% Setup Parmeters
clear all; close all;

f = 27.12e6;    %frequency of operation
w = 2*pi*f;
Rl = 500;       % Load on the secondary coil 
%(Rl ~~ for K<=0.2, <835.7 for K=0.3, <60.5 for K=1)
Ql = 20;        % Quality Factor for class E amplifier
K = 0.1;        % Coupling Coefficient
Ciss = 14.5e-12;% EPC2037: 14.5e-12; EPC8002: 20.5e-12
Coss = 17e-12;  % EPC2037: 17e-12; EPC8002: 13e-12
Crss = 1e-12;   % EPC2037: 1e-12; EPC8002: 1e-12
Co = Coss - Crss;% Typical Output Capacitance of the Switch
    
Lp = 507e-9;    % Primary coil Inductance
Rp = 1.46;      % Primary coil parasitic series Resistance
%Cpp = 2.8e-12; % Primary coil parasitic parallel Capacitance

Ls = 474e-9;    % Secondary coil Inductance
Rs = 1.13;      % Secondary coil parasitic series Resistance
%Cps = 2.5e-12; % Secondary coil parasitic parallel Capacitance
Cs = 75e-12;    % Capacitor for parallel resonance on the secondary side

Lc = 10e-6;     % Choke Inductance
Lcsi = 0.2e-9;  % Common Source Inductance
Rg = 1;         % Switch Gate Resistance
Lcable = 120e-9;% Cable Inductance
Rcable = 0.56;  % Cable Resistance

%% Transformer T - model
M = K*(Ls*Lp)^0.5;  % Mutual Inductance
L1 = Lp - M;        
L2 = Ls - M;
Z1 = Rs + j*w*L2 + Rl/(1 + j*w*Rl*Cs);
Z2 = Rp + j*w*L1 + j*w*M*Z1/(Z1 + j*w*M);
Z = Z2 + Rcable + j*w*Lcable % Impedence as seen from the primary side

%% Calculate Class E passive components parameters
C1 = 0.1836/(w*real(Z)) - Co
L = Ql*real(Z)/w - imag(Z)/w
C = 1/((Ql-1.1525)*w*real(Z))

%% Writing into the file
fileID = fopen('ClassE_param.txt','w');

% fprintf(fileID,'.tran 0 36.8u 27.6u\n');
% fprintf(fileID,';ac dec 4 1K 1G\n');
fprintf(fileID,'.four %dHz I(R1)\n\n',f);

fprintf(fileID,'.param f=%d Rl=%d Ql=%d\n', f,Rl,Ql);
fprintf(fileID,';.param C1=%d L=%d C=%d\n',C1,L,C);        % Class E Param
fprintf(fileID,'.param C1=1.42e-10 L=1.84e-07 C=4.8e-11\n');        % Class E Param
fprintf(fileID,'.param Lp=%d Rp=%d \n',Lp,Rp);  % Primary Coil Param
fprintf(fileID,'.param Ls=%d Rs=%d \n',Ls,Rs);  % Secondary Coil Param
fprintf(fileID,'.param Lc=%d Lcsi=%d Rg=%d Cs=%d\n',Lc,Lcsi,Rg,Cs);% Other Param
fprintf(fileID,'.param Lcab=%d Rcab=%d\n',Lcable,Rcable);   
fprintf(fileID,'K Lp Ls %d\n',K);
fprintf(fileID,';.step param Rl 100 1100 100\n\n');

fprintf(fileID,'.meas Pin AVG (I(V1)*V(Vdd))\n');
fprintf(fileID,'.meas Pgate AVG (I(V2)*V(Vg))\n');
fprintf(fileID,'.meas Ptx AVG (I(Lp)*V(Vtx))\n');
fprintf(fileID,'.meas Pout AVG (I(R1)*V(Vout))\n');

fprintf(fileID,'.meas PTE_PA param Ptx/Pin\n');
fprintf(fileID,'.meas PTE_coil param Pout/Ptx\n');
fprintf(fileID,'.meas PTE param Pout/Pin\n');
fprintf(fileID,'.meas PAE param Pout/(Pin+Pgate)\n');
fclose(fileID);
