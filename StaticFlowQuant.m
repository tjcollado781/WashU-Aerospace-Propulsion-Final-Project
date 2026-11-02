%Static flow quantities
A_7 = A_15 + A_x
I_7 = I_x + I_15
mfi_7 = 2*ht_7*(mdot_7/I_7)^2
M_7 = MbxMfi(mfi_7)

T_7 = Tt_7*TrixM(M_7,gamma)
V_7 = M_7*sqrt(gamma*R_7*T_7)
rho_7 = mdot_7/(V_7*A_7);            % slug/ft^3
P_7 = rho_7*R_7*T_7
Pt_7 = P_7 * (1 + (gamma - 1)/2 * M_7^2)^(gamma/(gamma - 1));  % psf
