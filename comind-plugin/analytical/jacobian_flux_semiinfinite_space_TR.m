function J = jacobian_flux_semiinfinite_space_TR(r_1, r_2, r_3, omega, mua, c, mus, g)
% This is an implementation of the Jacobian
%for the semi-infinite space as in REF[1].
%
% REFS:
%   [1] Arridge, S. R. (1995). Photon-measurement density functions. 
%       Part I: Analytical forms. Applied Optics, 34(31), 7395-7409.


b = 1;
z_1 = ( ( 1 - g ) * mus ) ^ (-1);
gamma = mua * c;
kappa = c / ( mua + ( 1 - g ) * mus);

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

% Compute Jacobian for absorption REF[1] Eq. 27
SI_A = compute_SI(rho_12_neg, rho_23, sigma, kappa);
SI_B = compute_SI(rho_12_pos, rho_23, sigma, kappa);

J_alpha = 2 * z_23 / rho_23 * ( SI_A -  SI_B );

% Compute Jacobian for diffusion REF[1] Eq. 28
TI_A = compute_TI(rho_12_neg, rho_23, sigma, kappa);
W_A = compute_WI(rho_12_neg, rho_23, sigma, kappa);
TI_B = compute_TI(rho_12_pos, rho_23, sigma, kappa);
W_B = compute_WI(rho_12_pos, rho_23, sigma, kappa);

J_nu =   2 * ( - z_12_neg * TI_A ...
             + z_23 * ( rho_tilde_12_neg * rho_tilde_23' ) * W_A )...
       - 2 * ( - z_12_pos * TI_B ... 
             + z_23 * ( rho_tilde_12_pos *  rho_tilde_23' ) *  W_B );

J = [J_alpha, J_nu];

end

function W = compute_WI(x, y, sigma, kappa)
% Implements REF[1] Eq. Table 2, W^{(\Gamma)}(x,y) for Time-Domain

W = 1 / ( 4 * ( 2 * pi ) ^ 3 * kappa ) ...
    * ( 3 + 3 * sigma * y + sigma ^ 2 * y ^ 2 ) / y ^ 4 ...
    * ( 1 + sigma * x ) / x ^ 2 ...
    * exp( - sigma * ( x + y ) );


end

function T = compute_TI(x, y, sigma, kappa)
% Implements REF[1] Eq. Table 2, T^{(\Gamma)}(x,y) for Time-Domain

T = 1 / ( 4 * ( 2 * pi ) ^ 3 * kappa ) ...
    * ( 1 + sigma * y ) / y ^ 3 ...
    * ( 1 + sigma * x ) / x ^ 3 ...
    * exp( - sigma * ( x + y ) );

end

function S = compute_SI(x, y, sigma, kappa)
% Implements REF[1] Eq. Table 2, S^{(\Gamma)}(x,y) for Time-Domain
S = 1 / ( 4 * ( 2 * pi ) ^ 3 * kappa ) * ( 1 + sigma * y ) / ( x * y ^ 2 ) ...
    * exp( - sigma * ( x + y ) ); 

end

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
