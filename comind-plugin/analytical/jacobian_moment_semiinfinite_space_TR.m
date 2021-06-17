function [J] = jacobian_moment_semiinfinite_space_TR(r_1, r_2, r_3, omega, mua, c, mus, g, isMellin)
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

% correct z_1 position, assuming z_1 is the boundary of the tissue we place
% the origin of the source deeper in the tissue REF[3]
b = 1;
z_1 = ( ( 1 - g ) * mus ) ^ (-1);
gamma = mua * c;
kappa = c / ( 3 * ( mua + ( 1 - g ) * mus ) );

if ( r_1(3) ~= 0 ) || ( r_3(3) ~= 0)
    error(['This function only supports z = 0  as the position of the boundary and origin of coordinates.',...
           'Sources and detectors should be placed at the boundary'])
end


if ~isMellin
    error("Not Implemented")
else
    J = jacobian_mellin_moment_semiinfinite_space_TR(r_1, r_2, r_3, z_1, omega, gamma, kappa);
end

end

function [J] = jacobian_mellin_moment_semiinfinite_space_TR(r_1, r_2, r_3, z_1, omega, gamma, kappa)

rho_23 = norm(r_3 - r_2);
rho_12_neg = compute_rho_ab(r_1 - [0 0 z_1], r_2);
rho_12_pos = compute_rho_ab(r_1 + [0 0 z_1], r_2);

rho_tilde_12_neg = compute_rho_ab_tilde(r_1 - [0 0 z_1], r_2);
rho_tilde_12_pos = compute_rho_ab_tilde(r_1 + [0 0 z_1], r_2);
rho_tilde_23 = compute_rho_ab_tilde(r_2, r_3);

% TODO are these below supposed to be distances over z?
z_23 = abs(r_3(3) - r_2(3));
z_12_neg = r_2(3) - z_1;
z_12_pos = r_2(3) + z_1;

sigma = compute_sigma(gamma, omega, kappa);


% Compute Jacobian for absorption REF[1] Eq. 47 
ST_A = compute_ST(rho_12_neg, rho_23, kappa, sigma);
ST_B = compute_ST(rho_12_pos, rho_23, kappa, sigma);
J_alpha =  2 * z_23 / rho_23 * (ST_A - ST_B);

% Compute Jacobian for diffusion REF[1] Eq. 52
% taking into account that to obtain the first moment we need to do
% omega = 0
TT_A = compute_TT(rho_12_neg, rho_23, sigma, kappa);
WT_A = compute_WT(rho_12_neg, rho_23, sigma, kappa);
TT_B = compute_TT(rho_12_pos, rho_23, sigma, kappa);
WT_B = compute_WT(rho_12_pos, rho_23, sigma, kappa);

J_nu =     - z_12_neg * TT_A + z_23 * ( rho_tilde_12_neg * rho_tilde_23' ) * WT_A ...
       - ( - z_12_pos * TT_B + z_23 * ( rho_tilde_12_pos * rho_tilde_23' ) * WT_B );
       

% Return as TODO
J = [J_alpha, J_nu];


end

function [J] = jacobian_temporal_moment_semiinfinite_space_TR(r_1, r_2, r_3, omega, gamma, kappa)

rho_13 = norm(r_3 - r_1);
rho_23 = norm(r_3 - r_2);
z_23 = rho_23(3);

rho_12_neg = compute_rho_ab_neg(r_1, r_2);
rho_12_pos = compute_rho_ab_pos(r_1, r_2);

sigma = compute_sigma(gamma, omega, kappa);

% Compute Jacobian for absorption REF[1] Eq. 48 
J_alpha =  z_23 / rho_23 * 1 / compute_G(r_3, r_1, sigma) ...
           * ( compute_ST(rho_12_neg, rho_23, kappa, sigma) - compute_ST(rho_12_pos, rho_23, kappa, sigma) ... 
               - compute_t(rho_13) * ( compute_SI(rho_12_neg, rho_23, kappa, sigma) - compute_ST(rho_12_pos, rho_23, kappa, sigma) ) );

% Compute Jacobian for diffusion REF[1] Eq. 52
J_nu = 1 / compute_G(r_3, r_1, sigma)...
       * ( jacobian_mellin_moment_semiinfinite_space_TR(z_23, rho_23, rho_12_neg, rho_12_pos)...
           - compute_t(rho_13) * jacobian_intensity_semiinfinite_space_TR(r_1, r_2, r_3) );
        

% Return as TODO
J = [J_alpha, J_nu];
end

function S = compute_ST(x, y, kappa, sigma)
% Implements REF[1] Eq. 43
    S = 1 / ( 8 * ( 2 * pi ) ^ 3 * kappa ^ 2 * sigma )...
        * ( x + sigma * ( x + y ) * y ) / ( x * y ^ 2 )...
        * exp( - sigma * ( x + y ) );
end

function T = compute_TT(x, y, sigma, kappa)
% Implements REF[1] Eq. 49.
% We assume that the hat comes from the fact that we can use non-stationary
% or stationary moments. We are interested in the stationary moments, i.e.
% the moments computed on the full TPSF, i.e. omega = 0. 
% Eq. 49 defines T but not T_hat. In addition Eq. 49 has a typo, inside 
% exp() there should be a sigma_0 not sigma.
T = 1 / ( 8 * ( 2 * pi ) ^ 3 * kappa ^ 2 )...
    * ( ( x ^ 2 + y ^ 2 ) + sigma * ( x ^ 2 * y + x * y ^ 2 ) ) / ( x * y ) ^ 3 ...
    * exp( - sigma * ( x + y ) );

end

function W = compute_WT(x, y, sigma, kappa)
% Implements REF[1] Eq. 50.
% We assume that the hat comes from the fact that we can use non-stationary
% or stationary moments. We are interested in the stationary moments, i.e.
% the moments computed on the full TPSF. Eq. 50 defines W but not W_hat.
% In addition Eq. 50 has a typo, inside exp() there should be a sigma_0 not
% sigma.

W = 1 / ( 8 * ( 2 * pi ) ^ 3 * kappa ^ 2 )...
    * ( 3 * ( x ^ 2 + y ^ 2 ) + sigma * ( x * y ^ 2 + 3 * x ^ 2 * y + y ^ 3 ) ...
        + sigma ^ 2 * ( x ^ 2 * y ^ 2 + x * y ^ 3 ) )...
    / ( x ^ 2 * y ^ 4 ) * exp( - sigma * ( x + y ) ); 

end

function t = compute_t(rho, kappa, sigma)
% As per REF[1] Eq. 44 and REF[2] Eqs. 4.1 and 4.2 in table 4
% we realise that REF[1] Eq. 44 is also valid for the semi-infinite space
t = 1 / ( 2 * kappa ) * rho ^ 2 / ( 1 + sigma * rho );
end

function G = compute_G(r_a, r_b, sigma)
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

rho_ab = r_b - r_a;

G = rho_ab(3) * ( 1 + sigma * norm(rho_ab) ) * exp( - sigma * norm(rho_ab) )...
    / ( ( 2 * pi ) ^ ( 3 / 2 ) * norm(rho_ab) ^ 3 );

end

function S = compute_SI(x, y, kappa, sigma)
% Implements REF[1] Eq. Table 2, S^{(\Gamma)}(x,y) for Frequency-Domain

S = 1 / ( 4 * ( 2 * pi ) ^ 3  * kappa ) * ( 1 + sigma * y ) / ( x * y ^ 2 )...
    * exp( - sigma * ( x + y ) );

end

% - Vector from position a to b: rho_vec_ab = r_b - r_a
% - Distance from a to b: rho_ab = norm(v_rho_ab)
% - Unit vector from position a to b: rho_tilde_ab = rho_vec_ab / rho_ab 
% - z_1 = b * ( ( 1 - f_hat ) * mu_s ) ^ (-1) from REF[3]

function rho_tilde = compute_rho_ab_tilde(r_a, r_b)
    rho_tilde = ( r_b - r_a ) / norm( r_b - r_a );
end

function sigma = compute_sigma(gamma, omega, kappa)
    % Implements REF[1] Table 1 row 1
    sigma = ((gamma + 1i*omega) / kappa) ^ (1/2);
end

function h = compute_h_ab(r_a, r_b)
    % Implements REF[1] Table 1 row 2 and 3
    % h = ( ( x_b - x_a ) ^ 2 + ( y_b - y_a ) ^ 2 ) ^ ( 1 / 2 );
    h = ( ( r_b(1) - r_a(1) ) ^ 2 + ( r_b(2) - r_a(2) ) ^ 2 ) ^ ( 1 / 2 );
end

function rho = compute_rho_ab(r_a, r_b)
    % Implements REF[1] Table 1 row 2 and 3
    % rho_ab_neg = (h_ab ^ 2 + ( z_b - z_a ) ^ 2 ) ^ ( 1 / 2 )
    h_ab = compute_h_ab(r_a, r_b);
    rho = (h_ab ^ 2 + ( r_b(3) + r_a(3) ) ^ 2 ) ^ ( 1 / 2 );
end
