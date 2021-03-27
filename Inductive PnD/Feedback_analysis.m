%% Instructions
% This file outputs the efficiency values of Class-E amplifier based on the corrected version of M. Baker et. al.
% Please go through the paper first.
% Both ideal and theoretical values are measured. All the parameters are programmable. 
% Make sure all the impedence values are legible.

%% Setup Parmeters
clear all; close all;

f = 27.12e6;    %frequency of operation
w = 2*pi*f;
Rac = 500;      % Load on the secondary coil 
%(Rac ~~ for K<=0.2, <835.7 for K=0.3, <60.5 for K=1)
%Ql = 50;        % Loaded quality Factor for class E amplifier
K = 0.12;       % Coupling Coefficient
Ciss = 14.5e-12;% EPC2037: 14.5e-12; EPC8002: 20.5e-12
Coss = 17e-12;  % EPC2037: 17e-12; EPC8002: 13e-12
Crss = 1e-12;   % EPC2037: 1e-12; EPC8002: 1e-12
Co = Coss-Crss; % Typical Output Capacitance of the Switch

Lp = 507e-9;    % Primary coil Inductance
Rp = 1.46;      % Primary coil parasitic series Resistance
%Cpp = 2.8e-12; % Primary coil parasitic parallel Capacitance

Ls = 474e-9;    % Secondary coil Inductance
Rs = 1.13;      % Secondary coil parasitic series Resistance
%Cps = 2.5e-12; % Secondary coil parasitic parallel Capacitance
Cs = 1/(w*w*Ls);    % Capacitor for parallel resonance on the secondary side
% 68e-12 for RB751SM-40; 56 for PMEG2010AEB

Lc = 10e-6;     % Choke Inductance
Lcsi = 0.2e-9;  % Common Source Inductance
Rg = 1;         % Switch Gate Resistance
Lcable = 120e-9;% Cable Inductance
Rcable = 0.56;  % Cable Resistance
Rl = 2*Rac;

%% Feedback modelling
Cp = 1/(w*w*Lp);     
M = 0.12*(Ls*Lp)^0.5;  % Mutual Inductance

Z1 = Rp + j*w*Lp;
Z2 = Rs + j*w*Ls + Rac/(1+j*w*Rac*Cs);

L = -w*w*M*M/(Z1*Z2);   
Zin = Z1*(1-L);
Z = Zin + Rcable + j*w*Lcable

Q1 = w*Lp/Rp;
Q2 = w*Ls/Rs;
Ql = w*Rac*Cs;
Q2l = Q2*Ql/(Q2+Ql);
Kc = (Q1*Q2l)^-0.5;

%% Plots
figure();
xlabel('coupling coefficient (K)');
ylabel('PTE');
hold on;
for K = 0.12
  Ql_opt = Q2/((1+K*K*Q1*Q2)^0.5);
  Rac_opt = Ql_opt/(w*Cs);
  
  PTE_coil = K*K*Q1*Q2l/(1+K*K*Q1*Q2l);
  PTE_sec = Q2/(Q2 + Ql);
  PTE = PTE_coil*PTE_sec;
  PTE_max = K*K*Q1*Q2/((1+(1+K*K*Q1*Q2)^0.5)^2);
  plot(K,100*PTE_max,'bo');
  plot(K,100*PTE,'ro');
end
filename='Data.csv';
PTE_sim=csvread(filename,1,1,[1 0 100 0]);
plot(0.01:0.01:1,PTE_sim,'go');

PTE_sim = zeros(100,1);

K = 0.12;
figure();
xlabel('Load resistance (Ohms)');
ylabel('PTE');
hold on;
% for Rac = 100:10:5000
%   Q1 = w*Lp/Rp;
%   Q2 = w*Ls/Rs;
%   Ql = w*Rac*Cs;
%   Q2l = Q2*Ql/(Q2+Ql);
%   Ql_opt = Q2/((1+K*K*Q1*Q2)^0.5);
%   Rac_opt = Ql_opt/(w*Cs);
%   
%   PTE_coil = K*K*Q1*Q2l/(1+K*K*Q1*Q2l);
%   PTE_sec = Q2/(Q2 + Ql);
%   PTE = PTE_coil*PTE_sec;
%   PTE_max = K*K*Q1*Q2/((1+(1+K*K*Q1*Q2)^0.5)^2);
%   plot(Rac,100*PTE_max,'bo');
%   plot(Rac,100*PTE,'ro');
% end