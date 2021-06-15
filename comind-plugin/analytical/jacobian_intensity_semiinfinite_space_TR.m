function J = jacobian_intensity_semiinfinite_space_TR(TODO)
% This is an implementation of the Jacobian
%for the semi-infinite space as in REF[1].
%
% REFS:
%   [1] Arridge, S. R. (1995). Photon-measurement density functions. 
%       Part I: Analytical forms. Applied Optics, 34(31), 7395-7409.

% Compute Jacobian for absorption REF[1] Eq. 27
J_alpha = 2 * z_23 / rho_23...
          * ( compute_SI(rho_12_neg, rho_23, omega) -  compute_SI(rho_12_pos, rho_23, omega) );

% Compute Jacobian for diffusion REF[1] Eq. 28
J_nu = 2 * ( - z_12_neg * compute_TI(rho_12_neg, rho_23, omega)...
             + z_23 * (rho_tilde_12_neg * rho_tilde_23) * compute_W(rho_12_neg, rho_23, omega) )...
       - 2 * ( - z_12_pos * compute_TI(rho_12_pos, rho_23, omega) ... 
               + z_23 * ( rho_tilde_12_neg *  rho_tilde_23 ) * compute_W(rho_12_pos, rho_23, omega) );

J = [J_alpha, J_nu];

end

function W = compute_W(x, y, t)
% Implements REF[1] Eq. Table 2, W^{(\Gamma)}(x,y) for Time-Domain

W = ( 3 * x / ( y ^ 4 * t )...
      + ( x + y ) ^ 2 * ( 3 * x ^ 2 * x * y + y ^ 2 ) / ( x ^ 2 + y ^ 3 )...
        * 1 / ( 2 * kappa * t ^ 2 )...
      + ( x + y ) ^ 4 / ( x * y ^ 2 ) * 1 / ( 4 * kappa ^ 2 * t ^ 3 ))...
    * 1 / ( 4 * kappa * ( 2 * pi ) ^ 3 ) * compute_f(x + y, t);


end

function T = compute_TI()
% Implements REF[1] Eq. Table 2, T^{(\Gamma)}(x,y) for Time-Domain

T = ( ( 1 / x + 1 / y ) ^ 3 * ( x * y ) / ( 2 * kappa * t ^ 2 )...
      + ( 1 / x ^ 3  + 1 / y ^ 3 ) * 1 / t )...
    * 1 / ( 4 * kappa * ( 2 * pi ) ^ 3 ) * compute_f(x + y, t);

end

function S = compute_SI(x, y, t, kappa, gamma)
% Implements REF[1] Eq. Table 2, S^{(\Gamma)}(x,y) for Time-Domain

S = 1 / ( 2 * ( 2 * pi ) ^ 3 )...
    * ( 1 / y ^ 2 + ( x + y ) / ( 2 * kappa * t ) * ( 1 / x + 1 / y ) )...
    * compute_f(x + y, t, gamma);

end

function f = compute_f(x, t, kappa, gamma)
% Implements REF[1] Eq. Table 2, f(x) for Time-Domain
 f = ( 2 * kappa * t ) ^ (- 3 / 2)...
     * exp( - gamma * t - x ^ 2 / ( 4 * kappa * t) );

end
