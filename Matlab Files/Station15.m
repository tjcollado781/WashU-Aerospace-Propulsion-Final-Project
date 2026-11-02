%% PROPULSION PROJECT PART 1

clc 
clear
close all

% parameters to set // iterate through
A0 = 4;
phiCh = 2;

% constants
gamma = 1.40;
Mwt_O2 = 32.00;
Mwt_N2 = 28.16;
Mwt_H2 = 2.00;
delta_hf = 1.30E9;
a = 8.80E-5;
b = 540;
Mach4Recovery = 0.669;
g0 = 32.17;

% Set up dummy freestream values %

Pt_Freestream = 22738; % DUMMY
Tt_Freestream = 1638;  % DUMMY

% solving for station 15 %

% continuity
Pt_inlet = Pt_Freestream*Mach4Recovery;
Tt_inlet = Tt_Freestream;
R_inlet = 1716;       % DUMMY
Cp_inlet = 6007;      % DUMMY
phi_inlet = 0;        % DUMMY
Beta_inlet = 3.762;   % DUMMY
md_inlet_H2 = 0.000;  % DUMMY
md_inlet_O2 = 0.8038; % DUMMY
md_inlet_N2 = 2.6615; % DUMMY

fracs_inlet = MassFracs(phi_inlet,Beta_inlet);
Mbp_inlet = MWT_yHyOyN(fracs_inlet(1),fracs_inlet(2),fracs_inlet(3)); % Molecular Weight Total

% givens for 15
M15_Mach4 = 0.2;

% pressures
Pt15 = Pt_inlet;
P15 = Pt15*PrixM(M15_Mach4,gamma);

% temperatures
Tt15 = Tt_inlet;
T15 = Tt15*TrixM(M15_Mach4,gamma);

% mass flows + fracs
phi15 = phi_inlet;
Beta15 = Beta_inlet;

Mbp15 = Mbp_inlet;

md15_H2 = md_inlet_H2;
md15_O2 = md_inlet_O2;
md15_N2 = md_inlet_N2;

md15 = md15_H2 + md15_O2 + md15_N2;

% gas constants

R15 = R_inlet;
Cp15 = Cp_inlet;

% velocity
V15 = M15_Mach4*sqrt(gamma*R15*Tt15);

% area
Rho15 = P15/(R15*T15);
A15 = md15/(V15*Rho15);

% impulse
I = md15*V15 + P15*A15;

% enthalpy
ht_15 = Cp15*Tt15;


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