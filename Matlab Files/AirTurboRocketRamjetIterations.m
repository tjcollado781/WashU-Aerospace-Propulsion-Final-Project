function outputs = AirTurboRocketRamjetIterations(alt,PI_f,phi) 
%Aerospace Propulsion - Semester Project
%Part 1 - Cruise Case (Mach 4.0, Altitude of 60000 ft)
%Compute m0 = m15 from A0.

%Enviromental Conditions (change between interations)
if alt == 60000 
%Cruise Case
    M_0 = 4.0;  %[mach]
    ThrustTarget = 5000; %[lbf]
    Cd_cowl = 0.4
    inletRecovery = 0.669;
    M_15 = 0.2;  %[mach]
else
%Launch Case
    M_0 = 0.8;  %[mach]
    ThrustTarget = 4000; %[lbf]
    Cd_cowl = 0.2
    inletRecovery = 0.96;
    M_15 = 0.2;  %[mach]
end

%Constants - (remain same throughout all iterations)
gamma = 1.4;
Mwt_O2 = 32.00;
Mwt_N2 = 28.16;
Mwt_H2 = 2.00;
R_universal = 49710; %[ft*lbf/kmol*(slug/kg)*R]
delta_hf = 1.30E9;
a = 8.80E-5;
b = 540;
g0 = 32.17;
Cfg = 0.96;


areaStep = 0.01;
A_0 = areaStep;


%Calculations from Freestream to Inlet
atmProperties = usAtm(alt);
P_0 = atmProperties(2);                    %[lbf/ft^2]
T_0 = atmProperties(1) + 459.67;                     %[R]
rho_0 = atmProperties(3);              %[slug/ft^3]
V_0 = M_0*atmProperties(4);                       %[ft/s]

Tt_0 = T_0*(1+((gamma-1)*M_0^2)/2);
Pt_0 = P_0*(Tt_0/T_0)^(gamma/(gamma-1));
q_0 = 0.5*rho_0*V_0^2;             %[lbf]

while true

A_cowl = 0.1*A_0;                  %[ft]
mDot_0 = rho_0*V_0*A_0;
D_cowl = q_0*Cd_cowl*A_cowl;
D_ram = mDot_0*V_0;


phi_1 = 0;        
Beta_1 = 3.762;  

Pt_1 = Pt_0*inletRecovery;
Tt_1 = Tt_0;
fracs_1 = MassFracs(phi_1,Beta_1);
Mbp_1 = MWT_yHyOyN(fracs_1(1),fracs_1(2),fracs_1(3)); % Molecular Weight Total
R_1 = R_universal/Mbp_1
Cp_1 = R_1*(gamma/(gamma-1))
md_1_H2 = fracs_1(1)*mDot_0;  
md_1_O2 = fracs_1(2)*mDot_0; 
md_1_N2 = fracs_1(3)*mDot_0; 
md_1 = md_1_H2 + md_1_O2 + md_1_N2;

% Fan
eta_f = 0.95;
R_fan = R_1;
Cp_fan = Cp_1;

Pt_fan = PI_f * Pt_1;
Tt_fan = Tt_1*(1+(1/eta_f)*(PI_f^((gamma-1)/gamma)-1));
W_fan = md_1*Cp_fan*(Tt_fan-Tt_1)

%Station 15
% pressures
Pt_15 = Pt_fan;
P_15 = Pt_15*PrixM(M_15,gamma);

% temperatures
Tt_15 = Tt_fan;
T_15 = Tt_15*TrixM(M_15,gamma);

% mass flows + fracs
phi_15 = phi_1;
Beta_15 = Beta_1;

Mbp_15 = Mbp_1;

md_15_H2 = md_1_H2;
md_15_O2 = md_1_O2;
md_15_N2 = md_1_N2;

md15 = md_15_H2 + md_15_O2 + md_15_N2;

% gas constants

R_15 = R_fan;
Cp_15 = Cp_fan;

% velocity
V_15 = M_15*sqrt(gamma*R_15*Tt_15);

% area
Rho_15 = P_15/(R_15*T_15);
A_15 = md15/(V_15*Rho_15);

% impulse
I_15 = md15*V_15 + P_15*A_15;

% enthalpy
ht_15 = md15*Cp_15*Tt_15;

%Rocket Motor 
Tt_ch = 540;
Pt_ch = 200000;
Cp_H2_ch = 87010;
Cp_O2_ch = 5438;
beta_ch = 0;
m_fuel = 1;
yONH_ch = MassFracs(phi,beta_ch);
a = 8.80*10^-5;
b = 540;
yH_ch = yONH_ch(1);
yO_ch = yONH_ch(2);
yN_ch = yONH_ch(3);

mDot_ch_H2 = yH_ch*m_fuel;
mDot_ch_O2 = yO_ch*m_fuel;
mDot_ch_N2 = yN_ch*m_fuel;

MWT_ch = MWT_yHyOyN(yH_ch,yO_ch,yN_ch)

R_ch = R_universal/MWT_ch
Cp_ch = R_ch*(gamma/(gamma-1))
ht_ch = (1/m_fuel)*(mDot_ch_H2*Cp_H2_ch*Tt_ch + mDot_ch_O2*Cp_O2_ch*Tt_ch)

Tt_b = Ttbrn_yHyOyNhi(yH_ch, yO_ch, yN_ch, ht_ch, gamma)
n_b = 1-a*(Tt_b-b)

mch_brnH2 = 69420;  % DONT ACC NEED THIS

% Throat, dont care
Tt_th = Tt_b
Pt_th = Pt_ch
M_th = 1
mDot_th = m_fuel
mft_th = MftxM(M_th)
A_th = (mDot_th*sqrt(R_ch*Tt_th))/(Pt_th*mft_th)

% turbine exit wahoo!


M_x = 1;

P_x = P_15;
Pt_x = (1/PrixM(M_x))*P_x;

eta_tau = 0.9;
eta_mu = 0.99;


PI_T = Pt_x/Pt_ch;
tau = 1-eta_tau*(1-PI_T^((gamma-1)/gamma));
Tt_x = tau*Tt_b;
SpecificTurbinePower = Cp_ch*(Tt_b-Tt_x)
md_x = W_fan/(eta_mu*SpecificTurbinePower)
T_x = TrixM(M_x)*Tt_x;

beta_x = beta_ch;
phi_x = phi;
yONH_x = MassFracs(phi_x,beta_x);

yH_x = yONH_x(1);
yO_x = yONH_x(2);
yN_x = yONH_x(3);

mDot_x_H2 = yH_x*md_x;
mDot_x_O2 = yO_x*md_x;
mDot_x_N2 = yN_x*md_x;

MWT_x = MWT_yHyOyN(yH_x,yO_x,yN_x);

R_x = R_universal/MWT_x;
Cp_x = R_x*(gamma/(gamma-1));
rho_x = (P_x*MWT_x)/(R_universal*T_x);
U_x = M_x*sqrt(gamma*R_x*T_x);
A_x = md_x/(rho_x*U_x);
I_x = md_x*U_x + (P_x-P_0)*A_x;
dh = ht_ch+(Cp_x*Tt_x+U_x^2/2)*md_x;



%Station 7
% Flow rates [sl/s]
mdot_H2_7 = mDot_x_H2 + md_1_H2;
mdot_O2_7 = mDot_x_O2 + md_1_O2;
mdot_N2_7 = mDot_x_N2 + md_1_N2;

% Feed total temperatures [R]
Tt_H2_7 = Tt_ch;
Tt_O2_7 = Tt_ch;
Tt_air_7 = Tt_1;

% Specific heats [ft^2/s^2/R]
Cp_H2_7 = Cp_H2_ch;
Cp_O2_7 = Cp_O2_ch;
Cp_air_7 = Cp_15;

%Total mass flow and mass fractions 
mdot_7 = mdot_H2_7 + mdot_O2_7 + mdot_N2_7;
y_H_7 = mdot_H2_7 / mdot_7;
y_O_7 = mdot_O2_7 / mdot_7;
y_N_7 = mdot_N2_7 / mdot_7;

%Mixture MW and gas constant R
Mb_7 = MWT_yHyOyN(y_H_7,y_O_7,y_N_7)
R_7 = R_universal / Mb_7;                        % ft^2/s^2/R
Cp_7 = R_7*(gamma/(gamma-1))       % ft^2/s^2/R

%Initial enthalpy h_i 
ht_7 = (1/mdot_7)*(mDot_x_H2*Cp_H2_7*Tt_H2_7 + mDot_x_O2*Cp_O2_7*Tt_O2_7+md_1*Cp_air_7*Tt_air_7) 


%Total temperature Tt7 after combustion 
Tt_7 = Ttbrn_yHyOyNhi(y_H_7,y_O_7,y_N_7, ht_7, gamma)

eta_b_7 = 1-((Cp_7*Tt_7 - ht_7) / (mdot_H2_7 * delta_hf));

%Static flow quantities
%rho_x = (P_x*MWT_x)/(R_universal*T_x);
% rho_7 = P_static_7 / (R_7 * T_static_7);            % slug/ft^3
% 
%        % ft/s
% A_7 = A_15+A_x;                 % ft^2
% M_7 = V_7 / sqrt(gamma * R_7 * T_static_7);
% Pt_7 = P_static_7 * (1 + (gamma - 1)/2 * M_7^2)^(gamma/(gamma - 1));  % psf
% T_7 = Tt_7*TrixM(M_7,gamma);
% P_7 = Pt_7*PrixM(M_7,gamma);
% V_7 = sqrt(2 * mdot_H2_7 * delta_hf * eta_b_7)/10; 

A_7 = A_15 + A_x
I_7 = I_x + I_15
mfi_7 = Cp_7*Tt_7*(mdot_7/I_7)^2
M_7 = MbxMfi(mfi_7)

T_7 = Tt_7*TrixM(M_7,gamma)
V_7 = M_7*sqrt(gamma*R_7*T_7)
rho_7 = mdot_7/(V_7*A_7);            % slug/ft^3
P_7 = rho_7*R_7*T_7
Pt_7 = P_7 * (1 + (gamma - 1)/2 * M_7^2)^(gamma/(gamma - 1));  % psf

%Enthalpy and combustion efficiency
Ht_7 = (Cp_7*Tt_7+V_7^2/2)*mdot_7
%I_7 = mdot_7*V_7 + (P_7-P_0)*A_7

m7_brnH2 = 69420;  % DONT ACC NEED THIS

%Nozzle
R_no = R_7
Cp_no = Cp_7
md_no = mdot_7
Tt_no = Tt_7
Pt_no = Pt_7
M_no = MxPri(P_0/Pt_no)
mft_no = MftxM(M_no)
A_no = (md_no*sqrt(R_no*Tt_no))/(Pt_no*mft_no)
T_no = Tt_no*TrixM(M_no,gamma);
P_no = Pt_no*PrixM(M_no,gamma);
V_no = M_no*sqrt(gamma*R_no*T_no)
F_no = md_no*V_no*Cfg



format longG
Freestream = [alt; M_0; P_0; Pt_0; T_0; Tt_0; q_0; Cd_cowl; A_0; V_0; mDot_0; A_cowl; D_cowl; D_ram]
Inlet = [inletRecovery; Pt_1; Tt_1; md_1; phi_1; Beta_1; md_1_H2; md_1_O2; md_1_N2; Mbp_1; R_1; Cp_1]
Station15 = [M_15;Pt_15;P_15;Tt_15;T_15;md15;A_15;V_15;I_15; ht_15;md_15_H2;md_15_O2;md_15_N2;Mbp_15;R_15;Cp_15]
Chamber = [Tt_ch; Pt_ch; m_fuel; phi; beta_ch; Cp_H2_ch; Cp_O2_ch; mDot_ch_H2; mDot_ch_O2; mDot_ch_N2; MWT_ch; R_ch; Cp_ch; ht_ch; n_b; Tt_th]
Throat = A_th
TurbineExit = [eta_tau; eta_mu; M_x; P_x; Pt_x; T_x; Tt_x; dh; md_x; phi_x; beta_x; mDot_x_H2; mDot_x_O2; mDot_x_N2; MWT_x; R_x; Cp_x; U_x; A_x; I_x]
Station7 = [mdot_7; mdot_H2_7; mdot_O2_7; mdot_N2_7; Mb_7; R_7; Cp_7; eta_b_7; Tt_7; Ht_7; I_7; A_7; rho_7; V_7; T_7; M_7; P_7; Pt_7]
Nozzle = [Mb_7; R_no; Cp_no; Tt_no; T_no; Pt_no;P_no;M_no;A_no;V_no;md_no;Cfg;F_no]

Tt_0 = T_0*TrixM(M_0)^-1;

Pt_0 = P_0*PrixM(M_0)^-1;

q_0 = 0.5*rho_0*V_0^2;             %[lbf]

A_cowl = 0.1*A_0;                  %[ft]

mDot_0 = rho_0*V_0*A_0;

D_cowl = q_0*Cd_cowl*A_cowl;

D_ram = mDot_0*V_0;


% Installed Thrust

F_installed = F_no - D_ram - D_cowl;
if round(F_installed) >= ThrustTarget
    break
end
    A_0 = A_0 + areaStep;
end



Cfg_no = (mdot_7 * V_7)/F_installed;

SFCeng = md_x*(32.174*3600)/F_installed

FNWA = F_installed/md_1 %Specific Thrust

ISP = 3600/SFCeng

FN_AMX = ThrustTarget/A_0

outputs = [SFCeng ISP FNWA A_no/A_0];

function resultsTable = usAtm(alt) 
%Outputs T in F, P in psf, densitry in slug/ft^3, and speed of sound it
%ft/s)
atmoData = [ -5000	76.84	17.554	27.45
0	59	14.696	23.77
5000	41.17	12.228	20.48
10000	23.36	10.108	17.56
15000	5.55	8.297	14.96
20000	-12.26	6.759	12.67
25000	-30.05	5.461	10.66
30000	-47.83	4.373	8.91
35000	-65.61	3.468	7.38
40000	-69.7	2.73	5.87
45000	-69.7	2.149	4.62
50000	-69.7	1.692	3.64
60000	-69.7	1.049	2.26
70000	-67.42	0.651	1.39
80000	-61.98	0.406	0.86
90000	-56.54	0.255	0.56
100000	-51.1	0.162	0.33
150000	19.4	0.02	0.037
200000	-19.78	0.003	0.0053
250000	-88.77	0	0.00065];

t = interp1(atmoData(:,1),atmoData(:,2),alt);

p = interp1(atmoData(:,1),atmoData(:,3),alt);

p = p*144;

d = interp1(atmoData(:,1),atmoData(:,4),alt)*10^-4;

s = sqrt(1.4*1715*(t + 459.6));

resultsTable = [t p d s];
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

function [Tri] = TrixM(M, g)
% Computes T/Tt for the given M.
  if nargin < 2, g = 1.4; end
  Tri = 1 / (1 + (g - 1) / 2 * M ^ 2);
end

function [Pri] = PrixM(M, g)
% Computes P/Pt for the given M.
  if nargin < 2, g = 1.4; end
  Pri = (1 + (g - 1) / 2 * M ^ 2) ^ (g / (1 - g));
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

function [Mb] = MbxMfi(MFIin, g)
% Computes the subsonic Mach number for the given mass flow function, ht*(mdot/I)^2

  if nargin < 2, g = 1.4; end
  eps = 0.0000001;
  itmax = 1000;

  Ma = 0;
  Fa = MfixM(Ma, g);
  Mb = 1;
  Fb = MfixM(Mb, g);
    
  for iter = 1: itmax
    if mod(iter, 2) == 1
      Mg = 0.5 * (Ma + Mb);
    else
      Mg = Ma + (sqrt(MFIin) - sqrt(Fa)) * (Mb - Ma) / (sqrt(Fb) - sqrt(Fa));
    end
    Fg = MfixM(Mg, g);
    
    if abs(Fg - MFIin) <= eps
      break;
    elseif Fg > MFIin
      Mb = Mg;
      Fb = Fg;
    else
      Ma = Mg;
      Fa = Fg;
    end
  end

  Mb = Mg;
end

function [Mfi] = MfixM(M, g)
% Computes ht*(mdot/I)^2 for the given M.
  if nargin < 2, g = 1.4; end
  Mfi = g ^ 2 / (g - 1) * M ^ 2 * (1 + (g - 1) / 2 * M ^ 2) / (1 + g * M ^ 2) ^ 2;
end

end