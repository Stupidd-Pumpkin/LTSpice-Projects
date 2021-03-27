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
Rac = 500;      % Load on the secondary coil 
%(Rac ~~ for K<=0.2, <835.7 for K=0.3, <60.5 for K=1)
Rl = 2*Rac;
Ql = 50;        % Loaded quality Factor for class E amplifier
K = 0.12;       % Coupling Coefficient
Ciss = 14.5e-12;% EPC2037: 14.5e-12; EPC8002: 20.5e-12
Coss = 17e-12;  % EPC2037: 17e-12; EPC8002: 13e-12
Crss = 1e-12;   % EPC2037: 1e-12; EPC8002: 1e-12
Co = Coss-Crss; % Typical Output Capacitance of the Switch

Lp = 507e-9;    % Primary coil Inductance; Old Coil: 380e-9
Rp = 1.46;      % Primary coil parasitic series Resistance; Old Coil: 0.32
%Cpp = 2.8e-12; % Primary coil parasitic parallel Capacitance

Ls = 474e-9;    % Secondary coil Inductance; Old Coil: 240e-9
Rs = 1.13;      % Secondary coil parasitic series Resistance; Old Coil: 0.28
%Cps = 2.5e-12; % Secondary coil parasitic parallel Capacitance
Cs = 56e-12;    % Capacitor for parallel resonance on the secondary side
    % 68e-12 for RB751SM-40; 56 for PMEG2010AEB; (Old coils: 127.5e-12 for PMEG2010AEB)

Lc = 10e-6;     % Choke Inductance
Lcsi = 0.2e-9;  % Common Source Inductance
Rg = 1;         % Switch Gate Resistance
Lcable = 120e-9;% Cable Inductance
Rcable = 0.56;  % Cable Resistance

Rlcp = 0.4;     % Choke parasitic resistance
Rc1p = 0.1;     % Shunt capacitor parasitic resistance
Rcp = 0.1;      % Series capacitor parasitic resistance
Rlp = 0.8;      % Inductor parasitic resistance
Rds = 0.55      % Switch on resistance

%% Transformer T - model
M = K*(Ls*Lp)^0.5;  % Mutual Inductance
L1 = Lp - M;
L2 = Ls - M;
Z1 = Rs + j*w*L2 + Rac/(1 + j*w*Rac*Cs);
Z2 = j*w*M*Z1/(Z1 + j*w*M);
Z = Z2 + Rcable + Rp + j*w*L1 + j*w*Lcable

%% Calculate Class E parameters
C1 = 0.1836/(w*real(Z)) - Co
L = Ql*real(Z)/w - imag(Z)/w
C = 1/((Ql-1.1525)*w*real(Z))

Vin = 0.8;
Po = 0.5768*Vin^2/real(Z2);;

%%Plots
figure();
xlabel('coupling coefficient (K)');
ylabel('PTE');
hold on;
for K = 0.01:0.01:1
M = K*(Ls*Lp)^0.5;  % Mutual Inductance
L1 = Lp - M;
L2 = Ls - M;
Z1 = Rs + j*w*L2 + Rac/(1 + j*w*Rac*Cs);
Z2 = j*w*M*Z1/(Z1 + j*w*M);
R = abs(Z2);
Ploss = (8*Rlcp/((pi^2+4)*R)) + ((pi^2+28)*0.55/(2*R*(pi^2+4))) + ...
        (Rc1p*(pi^2-4)/(2*R*(pi^2+4))) + (Rlp + Rcp + Rcable + Rp)/R;
PTE_PA = 100/(1+Ploss);
plot(K,PTE_PA,'bo');
end


%% Writing into the file
fileID = fopen('ClassE_rect_param.txt','w');

% fprintf(fileID,'.tran 0 36.8u 27.6u\n');
% fprintf(fileID,';ac dec 4 1K 1G\n');
fprintf(fileID,'.four %dHz I(Ls)\n\n',f);

fprintf(fileID,'.param f=%d Rl=%d Ql=%d\n', f,Rl,Ql);
fprintf(fileID,';.param C1=%d L=%d C=%d\n',C1,L,C);        % Class E Param
%fprintf(fileID,'.param C1=7.5e-11 L=1.78e-6 C=1.48e-11\n'); 
    %optimized for K=0.12 and RB751SM-40
fprintf(fileID,'.param C1=91e-12 L=750e-9 C=27e-12\n'); 
    % optimized for K=0.12 and PMEG2010AEB
%fprintf(fileID,'.param C1=1.4e-10 L=5.2e-7 C=3.2e-11\n'); 
    % optimized for K=0.1 and PMEG2010AEB
fprintf(fileID,'.param Lp=%d Rp=%d\n',Lp,Rp);  % Primary Coil Param
fprintf(fileID,'.param Ls=%d Rs=%d\n',Ls,Rs);  % Secondary Coil Param
fprintf(fileID,'.param Lc=%d Lcsi=%d Rg=%d Cs=%d\n',Lc,Lcsi,Rg,Cs);% Other Param
fprintf(fileID,'.param Lcab=%d Rcab=%d\n',Lcable,Rcable);   
fprintf(fileID,'K Lp Ls %d\n',K);
fprintf(fileID,';.step param k 0.02 0.2 0.02\n\n');

fprintf(fileID,'.meas Pin AVG -(I(V1)*V(Vdd))\n');
fprintf(fileID,'.meas Pgate AVG -(I(V2)*V(Vg))\n');
fprintf(fileID,'.meas Ptx AVG (I(Lp)*V(Vtx))\n');
fprintf(fileID,'.meas Prx AVG ((I(Ls)+I(Cs))*V(Vrx))\n');
fprintf(fileID,'.meas Pout AVG (I(R1)*(V(Vout+)-V(Vout-)))\n\n');

fprintf(fileID,'.meas PTE_PA param 100*Ptx/Pin\n');
fprintf(fileID,'.meas PTE_coil param 100*Prx/Ptx\n');
fprintf(fileID,'.meas PTE_rect param 100*Pout/Prx\n');
fprintf(fileID,'.meas PTE param 100*Pout/Pin\n');
fprintf(fileID,'.meas PAE param 100*Pout/(Pin+Pgate)\n');
fclose(fileID);
