%Aerospace Propulsion - Semester Project
%Part 1 - Cruise Case (Mach 4.0, Altitude of 60000 ft)
clc
close all
clear all

phi = 10;
A_0 = 4;
gamma = 1.4;

%Freestream/Station 0
alt = 60000;                    %[ft]
M_0 = 4.0;                     %[mach]
Cd_cowl = 0.4;                  %[]
P_0 =  150;                    %[lbf/ft^2]
T_0 = 390;                     %[R]
rho_0 = 2.2371e-04              %[slug/ft^3]
V_0 = 3872;                       %[ft/s]
Tt_0 = T_0*(1+((gamma-1)*M_0^2)/2)
Pt_0 = P_0*(Tt_0/T_0)^(gamma/(gamma-1))
q_0 = 0.5*rho_0*V_0^2;             %[lbf]
A_cowl = 0.1*A_0                  %[ft]
mDot_0 = rho_0*V_0*A_0
D_cowl = q_0*Cd_cowl*A_cowl
D_ram = mDot_0*V_0
Mwt_O2 = 32.00;
Mwt_N2 = 28.16;
Mwt_H2 = 2.00;
delta_hf = 1.30E9;
a = 8.80E-5;
b = 540;
Mach4Recovery = 0.669;
g0 = 32.17;

% Set up dummy freestream values %

Pt_0 = 22738; % DUMMY
Tt_0 = 1638;  % DUMMY

% solving for station 15 %

% continuity
Pt_inlet = Pt_0*Mach4Recovery;
Tt_inlet = Tt_0;
R_0 = 1716;       % DUMMY
Cp_0 = 6007;      % DUMMY
Beta_0 = 3.762;   % DUMMY
md_0_H2 = 0.000;  % DUMMY
md_0_O2 = 0.8038; % DUMMY
md_0_N2 = 2.6615; % DUMMY

fracs_0 = MassFracs(phi,Beta_0);
Mbp_0 = MWT_yHyOyN(fracs_0(1),fracs_0(2),fracs_0(3)); % Molecular Weight Total

% givens for 15
M15_Mach4 = 0.2;

% pressures
Pt15 = Pt_inlet;
P15 = Pt15*PrixM(M15_Mach4,gamma);

% temperatures
Tt15 = Tt_inlet;
T15 = Tt15*TrixM(M15_Mach4,gamma);

% mass flows + fracs
phi15 = phi;
Beta15 = Beta_0;

Mbp15 = Mbp_0;

md15_H2 = md_0_H2;
md15_O2 = md_0_O2;
md15_N2 = md_0_N2;

md15 = md15_H2 + md15_O2 + md15_N2;

% gas constants

R15 = R_0;
Cp15 = Cp_0;

% velocity
V15 = M15_Mach4*sqrt(gamma*R15*Tt15);

% area
Rho15 = P15/(R15*T15);
A15 = md15/(V15*Rho15);

% impulse
I = md15*V15 + P15*A15;

% enthalpy
ht_15 = Cp15*Tt15;

function [Pri] = PrixM(M, g)
% Computes P/Pt for the given M.
  if nargin < 2, g = 1.4; end
  Pri = (1 + (g - 1) / 2 * M ^ 2) ^ (g / (1 - g));
end

function [ yOHN ] = MassFracs( Phi, Beta )
% Computes the mass fractions of oxygen, hydrogen, and nitrogen atoms
% in the mixture.
% Inputs: Phi  = mixture equivalence ratio
%         Beta = molar ratio of nitrogen to oxygen atoms in the mixture
% Output: Yh = mass fraction of hydrogen atoms
%
% Assumed molecular weights:
%         O  = 16
%         H  =  1
%         N  = 14.082
%             (adjusted to account for nitrogen, argon, etc. in air.)

% Set the number of moles of each component, assuming one mole of oxygen atoms.
  nO = 1;
  nN = Beta * nO;
  nH = 2 * Phi * nO;

% Calculate the resulting mixture mass and the mass fration of hydrogen atoms
  mMix = 16 * nO + 14.082 * nN + 1 * nH;
  yO   = (16 * nO) ./ mMix;
  yH   = (1 * nH) ./ mMix;
  yN   = (14.082 * nN) ./ mMix;
  
  yOHN = [ yH, yO, yN ];
end

function [ MWT ] = MWT_yHyOyN(yH, yO, yN)
% Computes the mixture average molecular weight assuming the hydrogen-oxygen reaction goes to completion
% Inputs: yH = mass fraction of hydrogen atoms
%         yO = mass fraction of oxygen atoms
%         yN = mass fraction of inert atoms
% 
% Assumed molecular weights: 
%         O2  = 32
%         H2  =  2
%         H2O = 18
%         N2  = 28.164
%              (adjusted to account for nitrogen, argon, etc. in air.)

% Compute the equivalence ratio and then, for one half mole of O2 reactant,
% find the number of moles of each resulting consitutuent including water,
% and unburned hydrogen, oxygen, and nitrogen (inerts).
  Phi = 8 * (yH ./ yO);
  
  nH2O = min(Phi, 1);
  nH2u = max(Phi - 1, 0);
  nO2u = 0.5 * max(1 - Phi, 0);
  nN2  = 0.5 * (yN / 28.164) ./ (yO / 32);

% Find the mixture weight and average molecular weight.
  mH2O = 18 * nH2O;
  mH2u = 2 * nH2u;
  mO2u = 32 * nO2u;
  mN2  = 28.164 * nN2;
  
  MWT = (mH2O + mH2u + mO2u + mN2) ./ (nH2O + nH2u + nO2u + nN2);

end

function [ Ttbrn ] = Ttbrn_yHyOyNhi(yH, yO, yN, hi, g)
% Approximates the gas total temperature after combution of hydrogen
% including a crude estimate to account for dissociation.
% Inputs: yH  = mass fraction of hydrogen atoms
%         yO  = mass fraction of oxygen atoms
%         yN  = mass fraction of inert atoms
%         hi  = specific enthalpy of reactants, not including heats of formation [ft^2/s^2]
% 
% 'Combustion of hydrogen approximated by:
%         hpr = heat of combustion = 1.3e9 [ft^2/s^2]
%         a   = combustion efficiency/dissociation sensitivity = 8.8e-5
%         b   = combustion efficiency/dissociation constant = 540

% Set constants
  hpr = 1300000000;
  a = 0.000088;
  b = 540;
 
% Ensure mass fractions add to one and get the combusted gas properties.
  Ytot = yH + yO + yN;
  yH = yH / Ytot;
  yO = yO / Ytot;
  yN = yN / Ytot;

  MWT = MWT_yHyOyN(yH, yO, yN);
  R = 49710. / MWT;
  Cp = R * g / (g - 1);

% Compute the mass fraction of burned hydrogen and the mixture temperature.
  Yh2b = min(yH, yO / 8);
  
  Ttbrn = (hi + Yh2b * hpr * (1 + a * b)) / (Cp + a * Yh2b * hpr);

end

function [Tri] = TrixM(M, g)
% Computes T/Tt for the given M.
  if nargin < 2, g = 1.4; end
  Tri = 1 / (1 + (g - 1) / 2 * M ^ 2);
end

function [Mft] = MftxM(M, g)
% Computes mdot(RTt)^.5/(PtA) for the given M.
  if nargin < 2, g = 1.4; end
  Mft = g ^ 0.5 * M * (1 + (g - 1) / 2 * M ^ 2) ^ (-(g + 1) / (2 * (g - 1)));
end

function [M] = MxPri(PR, g)
% Computes M for the given P/Pt.
  if nargin < 2, g = 1.4; end
  M = (2 / (g - 1) * (PR ^ ((1 - g) / g) - 1)) ^ 0.5;
end

function [Ari] = ArixM(M, g)
% Computes A/A* for the given M.
  if nargin < 2, g = 1.4; end
  Ari = (2 * (1 + (g - 1) / 2 * M ^ 2) / (g + 1)) ^ ((g + 1) / (2 * (g - 1))) / M;
end



