%Aerospace Propulsion - Semester Project
%Part 1 - Cruise Case (Mach 4.0, Altitude of 60000 ft)
%Compute m0 = m15 from A0.
clc
close all
clear all

phi = 2;
A_0 = 4;
gamma = 1.4;

%Freestream/Station 0
alt = 60000;                    %[ft]
M_0 = 4.0;                     %[mach]
Cd_cowl = 0.4;                  %[]
P_0 =  150;                    %[lbf/ft^2]
T_0 = 390;                     %[R]
%rho_sl = 0.002377;              %[slug/ft^3]
%rho_0 = rho_sl*9.4919*10^-2;   %[slug/ft^3]
rho_0 = 2.2371e-04              %[slug/ft^3]
V_0 = 3872;                       %[ft/s]
Tt_0 = T_0*(1+((gamma-1)*M_0^2)/2)
Pt_0 = P_0*(Tt_0/T_0)^(gamma/(gamma-1))
q_0 = 0.5*rho_0*V_0^2;             %[lbf]
A_cowl = 0.1*A_0                  %[ft]
mDot_0 = rho_0*V_0*A_0
D_cowl = q_0*Cd_cowl*A_cowl
D_ram = mDot_0*V_0


function [Pri] = PrixM(M, g)
% Computes P/Pt for the given M.
  if nargin < 2, g = 1.4; end
  Pri = (1 + (g - 1) / 2 * M ^ 2) ^ (g / (1 - g));
end

function [Tri] = TrixM(M, g)
% Computes T/Tt for the given M.
  if nargin < 2, g = 1.4; end
  Tri = 1 / (1 + (g - 1) / 2 * M ^ 2);
end



