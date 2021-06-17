function J = jacobian_flux_time_semiinfinite_space_TR(r_1, r_2, r_3, t, mua, c, mus, g)
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
musp = ( 1 - g ) * mus;
z_1 = 1 / musp;
gamma = mua * c;
kappa = c / ( 3 * ( mua + musp ) );

% Compute positions
r_1_neg = r_1 + [0 0 z_1]; % By convention in the z+ side
r_1_pos = r_1 - [0 0 z_1]; % By convention in the z- side (mirror source)
rho_23 = norm(r_3 - r_2);
rho_12_neg = compute_rho_ab(r_1_neg, r_2);
rho_12_pos = compute_rho_ab(r_1_pos, r_2);

rho_tilde_12_neg = compute_rho_ab_tilde(r_1_neg, r_2);
rho_tilde_12_pos = compute_rho_ab_tilde(r_1_pos, r_2);
rho_tilde_23 = compute_rho_ab_tilde(r_2, r_3);

% For a detector s(r_3) = [0, 0, -1], this are the z-dim of the equivalent
% r_ab = r_b - r_a vectors
z_23 = (r_3(3) - r_2(3));
z_12_neg = (r_2(3) - r_1_neg);
z_12_pos = (r_2(3) - r_1_pos);

% Compute Jacobian for absorption REF[1] Eq. 27
SI_A = compute_SI(rho_12_neg, rho_23, t, kappa, gamma);
SI_B = compute_SI(rho_12_pos, rho_23, t, kappa, gamma);
J_alpha = 2 * z_23 / rho_23 * ( SI_A -  SI_B );

% Compute Jacobian for diffusion REF[1] Eq. 28
TI_A = compute_TI(rho_12_neg, rho_23, t, kappa, gamma);
W_A = compute_WI(rho_12_neg, rho_23, t, kappa, gamma);
TI_B = compute_TI(rho_12_pos, rho_23, t, kappa, gamma);
W_B = compute_WI(rho_12_pos, rho_23, t, kappa, gamma);

J_nu =   2 * (...
                (...
                    - z_12_neg * TI_A ...
                    + z_23 * ( rho_tilde_12_neg * rho_tilde_23' ) * W_A...
                )...
              - (...
                    - z_12_pos * TI_B ... 
                    + z_23 * ( rho_tilde_12_pos *  rho_tilde_23' ) *  W_B...
                )...
             );

J = [J_alpha, J_nu];

end

function W = compute_WI(x, y, t, kappa, gamma)
% Implements REF[1] Eq. Table 2, W^{(\Gamma)}(x,y) for Time-Domain

W = (...
        3 * x / ( y ^ 4 * t )...
        + ( x + y ) ^ 2 * ( 3 * x ^ 2 - 2 * x * y + y ^ 2 ) / ( x ^ 2 * y ^ 3 )...
        * ( 1 / ( 2 * kappa * t ^ 2 ) )...
        + ( x + y ) ^ 4 / ( x * y ^ 2) * 1 / ( kappa * kappa ^ 2 * t ^ 3)...
    )...
    * 1 / ( 4 * kappa * ( 2 * pi ) ^ 3 ) * compute_f(x + y, t, kappa, gamma);

end

function T = compute_TI(x, y, t, kappa, gamma)
% Implements REF[1] Eq. Table 2, T^{(\Gamma)}(x,y) for Time-Domain
T = ( ( 1 / x  + 1 / y ) ^ 3 * ( x * y / ( 2 * kappa * t ^ 2) )...
      + ( 1 / x ^ 3  + 1 / y ^ 3) * 1 / t )...
    * 1 / ( 4 * kappa * ( 2 * pi ) ^ 3 ) * compute_f(x + y, t, kappa, gamma);

end

function S = compute_SI(x, y, t, kappa, gamma)
    S = 1 / ( 2 * ( 2 * pi ) ^ 3 )...
        * ( 1 / y ^ 2 + ( x + y ) / ( 2 * kappa * t ) * ( 1 / x + 1 / y ) )...
        * compute_f(x + y, t, kappa, gamma);

end

function f = compute_f(x, t, kappa, gamma)
    f = ( 2 * kappa * t ) ^ ( - 3 / 2 )...
        * exp( - gamma * t - x ^ 2 / ( 4 * kappa * t ) );
end

function rho_tilde = compute_rho_ab_tilde(r_a, r_b)
    rho_tilde = ( r_b - r_a ) / norm( r_b - r_a );
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
    rho = ( h_ab ^ 2 + ( r_b(3) - r_a(3) ) ^ 2 ) ^ ( 1 / 2 );
end
