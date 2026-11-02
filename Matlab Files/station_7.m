
function station7_simplified()
clc; clear;

% Flow rates [sl/s]
mdot_H2 = 0.034;
mdot_O2 = 0.941;
mdot_N2 = 2.661;

% Feed total temperatures [R]
Tt_H2 = 540;
Tt_O2 = 540;
Tt_air = 1638;

% Specific heats [ft^2/s^2/R]
Cp_H2 = 87010;
Cp_O2 = 5438;
Cp_air = 6339;

% Molecular weights [lbm/mol]
MW_H2 = 2.00;
MW_O2 = 32.00;
MW_N2 = 28.164;

% Constants
hpr = 1.3e9;           % heat of reaction [ft^2/s^2]
a = 8.8e-5; b = 540;
gamma = 1.4;
T_static = 3020;       % R
P_static = 14860;      % psf

%Total mass flow and mass fractions 
mdot_total = mdot_H2 + mdot_O2 + mdot_N2;
y_H = mdot_H2 / mdot_total;
y_O = mdot_O2 / mdot_total;
y_N = mdot_N2 / mdot_total;

%Mixture MW and gas constant R
Mb = 1 / (y_H / MW_H2 + y_O / MW_O2 + y_N / MW_N2);
R = 49710 / Mb;                        % ft^2/s^2/R
Cp = gamma * R / (gamma - 1);        % ft^2/s^2/R

%Initial enthalpy h_i 
H2 = mdot_H2 * Cp_H2 * Tt_H2
O2 =mdot_O2 * Cp_O2 * Tt_O2
N2 = mdot_N2 * Cp_air * Tt_air
h_i = (H2+O2+N2)/mdot_total

%Total temperature Tt7 after combustion 
numerator = h_i +0.0187 * hpr * (1 + (a/R) * b);
denominator = Cp + a *0.0187 * hpr;
Tt7 = numerator / denominator;


%Enthalpy and combustion efficiency
ht7 = Cp * Tt7;
eta_b = 1-((ht7 - h_i) / (mdot_H2 * hpr));

%Static flow quantities
rho = P_static / (R * T_static);            % slug/ft^3
U = sqrt(2 * mdot_H2 * hpr * eta_b)/10;        % ft/s
A = mdot_total / (rho * U);                 % ft^2
M = U / sqrt(gamma * R * T_static);
Pt = P_static * (1 + (gamma - 1)/2 * M^2)^(gamma/(gamma - 1));  % psf

%Gross Thrust 
F_gross = 18687 % from station 9
% Ram Drag
D_cowl = 268.3157;
% Cowl (Additive) Drag
D_ram = 1.3416e+04;

% Installed Thrust
F_installed = F_gross - D_ram - D_cowl;
Cfg = (mdot_total * U)/F_installed ;

%RESULTS 
fprintf('\n--- Station 7 Step-by-Step Results ---\n');
fprintf('mdot_total     = %.4f sl/s\n', mdot_total);
fprintf('Mb             = %.2f lbm/mol\n', Mb);
fprintf('R              = %.2f ft^2/s^2/R\n', R);
fprintf('Cp             = %.2f ft^2/s^2/R\n', Cp);
fprintf('h_i            = %.2e ft^2/s^2\n', h_i);
fprintf('Tt7            = %.0f R\n', Tt7);
fprintf('ht7            = %.2e ft^2/s^2\n', ht7);
fprintf('eta_b          = %.2f %%\n', 100 * eta_b);
fprintf('rho            = %.2e sl/ft^3\n', rho);
fprintf('U              = %.0f ft/s\n', U);
fprintf('A              = %.2e ft^2\n', A);
fprintf('Mach           = %.3f\n', M);
fprintf('Pt             = %.0f psf\n', Pt);
fprintf('F              = %.1f lbf\n', F_installed);
fprintf('Cfg            = %.3f\n', Cfg);
end
