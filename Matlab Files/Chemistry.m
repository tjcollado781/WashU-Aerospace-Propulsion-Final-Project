%Aerospace Propulsion - Semester Project
%Part 1 - Cruise Case (Mach 4.0, Altitude of 60000 ft)
%For the selected phi, (and b=0) compute mass flows of H and O
clc
close all
clear all

Tt_ch = 540;
Pt_ch = 200000;
Cp_H2_ch = 87010;
Cp_O2_ch = 5438;
phi_ch = 2;
A_0 = 4;
gamma = 1.4;
beta = 0;
m_fuel = 0.172;
R_ch = 49710;
yONH_ch = MassFracs(phi_ch,beta);
a = 8.80*10^-5;
b = 540;
hpr = 1.3*10^9;
yH_ch = yONH_ch(1);
yO_ch = yONH_ch(2);
yN_ch = yONH_ch(3);

mDot_ch_H2 = yH_ch*m_fuel;
mDot_ch_O2 = yO_ch*m_fuel;
mDot_ch_N2 = yN_ch*m_fuel;

MWT_ch = MWT_yHyOyN(yH_ch,yO_ch,yN_ch)

R_ch = R/MWT_ch
Cp_ch = R*(gamma/(gamma-1))
ht_ch = (1/m_fuel)*(mDot_ch_H2*Cp_H2_ch*Tt_ch + mDot_ch_O2*Cp_O2_ch*Tt_ch)

Tt_b = Ttbrn_yHyOyNhi(yH_ch, yO_ch, yN_ch, ht_ch, gamma)
n_ch = 1-a*(Tt_ch-b)

Tt_th = Tt_b
Pt_th = Pt_ch
M_th = 1
mDot_th = m_fuel
mft_th = MftxM(M_th)
A_th = (mDot_th*sqrt(R_ch*Tt_th))/(Pt_th*mft_th)

P15 = 14794
Px = P15

Mx = MxPri((Px/Pt_th))
Tx = TrixM(Mx)*Tt_th
Ux = Mx*sqrt(gamma*R_ch*Tx)
Ax = ArixM(Mx)*A_th
P = 150
I = mDot_th*Ux + (Px-P)*Ax
Htx = (ht_ch+Ux^2/2)*0.172

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
