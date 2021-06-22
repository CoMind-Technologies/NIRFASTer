function J = jacobian_flux_time_semiinfinite_space_TD(r_1, r_2, r_3, t, mua, c, mus, g) 
% This is an implementation of the Jacobian
%for the semi-infinite space as in REF[1].
%
% REFS:
%   [1] Arridge, S. R. (1995). Photon-measurement density functions. 
%       Part I: Analytical forms. Applied Optics, 34(31), 7395-7409.

if r_1(3) ~= 0 || r_3(3) ~= 0
    error("Source and detector should be applied at the boundary z=0.")
end

% Compute optical parameters
[z_1, gamma, kappa] = compute_optical_parameters(c, mua, mus, g);

% Compute positions
[rho_23, rho_12_neg, rho_12_pos, rho_tilde_12_neg, ...
 rho_tilde_12_pos, rho_tilde_23, z_23, z_12_neg, z_12_pos] = compute_positions(r_1, r_2, r_3, z_1);

% Compute Jacobian for absorption REF[1] Eq. 27
S_neg = compute_S_TD(rho_12_neg, rho_23, t, kappa, gamma); %<ok>
S_pos = compute_S_TD(rho_12_pos, rho_23, t, kappa, gamma); %<ok>
J_alpha = 2 * z_23 / rho_23 * ( S_neg - S_pos );           %<ok>

% Compute Jacobian for diffusion REF[1] Eq. 28
T_neg = compute_T_TD(rho_12_neg, rho_23, t, kappa, gamma); %<ok>
W_neg = compute_W_TD(rho_12_neg, rho_23, t, kappa, gamma); %<ok>
T_pos = compute_T_TD(rho_12_pos, rho_23, t, kappa, gamma); %<ok>
W_pos = compute_W_TD(rho_12_pos, rho_23, t, kappa, gamma); %<ok>

J_nu =   2 * (... %<ok>
                (...  %<ok>
                    - z_12_neg * T_neg ... %<ok>
                    + z_23 * ( rho_tilde_12_neg * rho_tilde_23' ) * W_neg... %<ok>
                )... %<ok>
              - (... %<ok>
                    - z_12_pos * T_pos ...  %<ok>
                    + z_23 * ( rho_tilde_12_pos *  rho_tilde_23' ) *  W_pos... %<ok>
                )...
             );  %<ok>

J = [J_alpha, J_nu];

end


