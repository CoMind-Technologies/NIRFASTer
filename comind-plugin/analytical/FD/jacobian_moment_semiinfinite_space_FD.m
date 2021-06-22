function [J] = jacobian_moment_semiinfinite_space_FD(r_1, r_2, r_3, mua, c, mus, g, order, isMellin)
%JACOBIAN_MOMENT_SEMIINFINITE_SPACE_TR This is an implementation of the Jacobian
%for the semi-infinite space as in REF[1].
%
% NOTATION: REF[1] Section A
%
% - Position vector r = [x, y, z]
%   - r_1 source
%   - r_2 point of Jacobian evaluation
%   - r_3 detector
% - Vector from position a to b: rho_vec_ab = r_b - r_a
% - Distance from a to b: rho_ab = norm(v_rho_ab)
% - Unit vector from position a to b: rho_tilde_ab = rho_vec_ab / rho_ab 
% - sigma = ((gamma + i*omega) / kappa)^ (1/2) for omega > 0
% - sigma_0 = sigma(omega = 0)
% - h_ab = ( ( x_b - x_a ) ^ 2 + ( y_b - y_a ) ^ 2 ) ^ ( 1 / 2 )
% - rho_12_neg = (h_12 ^ 2 + ( z_2 - z_1 ) ^ 2 ) ^ ( 1 / 2 )
% - rho_12_pos = (h_12 ^ 2 + ( z_2 + z_1 ) ^ 2 ) ^ ( 1 / 2 )
% - z_1 = b * ( ( 1 - f_hat ) * mu_s ) ^ (-1) from REF[3]
% - mua, linear absorption coefficient
% - mus, linear scattering coefficient
% - mus_prime, reduced scattering coefficient, (1 - g) * mus
% - g, mean cosine of the scattering angle or anisotropy factor (0 <= g <= 1)
% - D, diffusion coefficient: D = 3 * ( mua + (1 - h) * mu_s ) ^ {-1}
% - c, speed of light in the tissue
% - gamma = mua*c
% - kappa = c / (3 * ( mua + mus_prime ) )
%
% REFS:
%   [1] Arridge, S. R. (1995). Photon-measurement density functions. 
%       Part I: Analytical forms. Applied Optics, 34(31), 7395-7409.
%   [2] Arridge, S. R., Cope, M., & Delpy, D. T. (1992). The theoretical 
%       basis for the determination of optical pathlengths in tissue: 
%       temporal and frequency analysis. 
%       Physics in Medicine & Biology, 37(7), 1531.
%   [3] Patterson, M. S., Chance, B., & Wilson, B. C. (1989). 
%       Time resolved reflectance and transmittance for the noninvasive 
%       measurement of tissue optical properties. 
%       Applied optics, 28(12), 2331-2336.
%

if ( r_1(3) ~= 0 ) || ( r_3(3) ~= 0)
    error(['This function only supports z = 0  as the position of the boundary and origin of coordinates.',...
           'Sources and detectors should be placed at the boundary'])
end

% Compute sigma for omega = 0 since we are interested in the moments when
% t -> inf, i.e. the full TPSF has been captured which in Fourier Domain
% translates into w --> 0.

omega = 0;

if order == 1
    if ~isMellin
        J = jacobian_temporal_moment_semiinfinite_space_FD(r_1, r_2, r_3, 0, mua, c, mus, g);
    else
        J = jacobian_mellin_moment_semiinfinite_space_FD(r_1, r_2, r_3, 0, mua, c, mus, g);
    end
else
    error("Not implemented")
end

end

function [J] = jacobian_mellin_moment_semiinfinite_space_FD(r_1, r_2, r_3, omega, mua, c, mus, g)

% Compute optical parameters
[z_1, gamma, kappa] = compute_optical_parameters(c, mua, mus, g);

% Compute positions
[rho_23, rho_12_neg, rho_12_pos, rho_tilde_12_neg, ...
 rho_tilde_12_pos, rho_tilde_23, z_23, z_12_neg, z_12_pos] = compute_positions(r_1, r_2, r_3, z_1);

sigma = compute_sigma(gamma, omega, kappa);

% Compute Jacobian for absorption REF[1] Eq. 47 
ST_A = compute_S_FD(rho_12_neg, rho_23, kappa, sigma);
ST_B = compute_S_FD(rho_12_pos, rho_23, kappa, sigma);
J_alpha =  2 * z_23 / rho_23 * (ST_A - ST_B);

% Compute Jacobian for diffusion REF[1] Eq. 52
% taking into account that to obtain the first moment we need to do
% omega = 0
TT_A = compute_T_FD(rho_12_neg, rho_23, sigma, kappa);
WT_A = compute_W_FD(rho_12_neg, rho_23, sigma, kappa);
TT_B = compute_T_FD(rho_12_pos, rho_23, sigma, kappa);
WT_B = compute_W_FD(rho_12_pos, rho_23, sigma, kappa);

J_nu =     - z_12_neg * TT_A + z_23 * ( rho_tilde_12_neg * rho_tilde_23' ) * WT_A ...
       - ( - z_12_pos * TT_B + z_23 * ( rho_tilde_12_pos * rho_tilde_23' ) * WT_B );
       
J = [J_alpha, J_nu];

end

function [J] = jacobian_temporal_moment_semiinfinite_space_FD(r_1, r_2, r_3, omega, mua, c, mus, g)
% Compute G^{E}_{half}.
% We depart from REF[1] Eq. 26 which defines G(\omega)^{\Lambda}_{half}
% By the definition
%   G(t)^{E}_{half} = \int^{inf}_{0} G(t)^{\Lambda}_{half} dt 
% Taking into account that
%    G(t)^{\Lambda}_{half} = 0 for all t < 0
% Then
%    G^{E}_{half} = \int^{inf}_{-inf} G(t)^{\Lambda}_{half} dt
% And from REF[1] Eq. 36 in Fourier space
%   \int^{inf}_{-inf} G(t)^{\Lambda}_half dt = 
%       \sqrt{2\pi} G(\omega)^{\Lambda}_{half} |_(\omega=0)
% So finally,
%   G^{E}_{half} = \sqrt{2\pi} G(\omega)^{\Lambda}_{half} |_(\omega=0)

% Compute optical parameters
[z_1, gamma, kappa] = compute_optical_parameters(c, mua, mus, g);

% Compute positions
[rho_23, rho_12_neg, rho_12_pos, ~, ...
 ~, ~, z_23, ~, ~] = compute_positions(r_1, r_2, r_3, z_1);

sigma = compute_sigma(gamma, 0, kappa);

% Compute Jacobian for absorption REF[1] Eq. 47 
ST_A = compute_S_FD(rho_12_neg, rho_23, kappa, sigma);
ST_B = compute_S_FD(rho_12_pos, rho_23, kappa, sigma);
G_E  = compute_G_FD(r_3, r_2, sigma);
J_alpha =  2 * z_23 / rho_23 * 1 / G_E * (ST_A - ST_B);

% TODO and still TODO many questions here
% Compute Jacobian for diffusion REF[1] Eq. 52
% taking into account that to obtain the first moment we need to do
% omega = 0
JT = jacobian_mellin_moment_semiinfinite_space_FD(r_1, r_2, r_3, omega, mua, c, mus, g);
JT_nu = JT(2);
rho_13 = norm( r_3 - r_1 );
t_Gamma = compute_t_Gamma_FD(rho_13, sigma, kappa);
J_Gamma = jacobian_flux_semiinfinite_space_FD(r_1, r_2, r_3, mua, c, mus, g, omega);
J_nu_Gamma = J_Gamma(2);

J_nu = 1 / G_E * (  JT_nu - t_Gamma * J_nu_Gamma );
       
J = [J_alpha, J_nu];

end


