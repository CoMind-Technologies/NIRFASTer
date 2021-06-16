function [J] = jacobian_moment_semiinfinite_space_TR(chi_3,r1,r2)
%JACOBIAN_MOMENT_SEMIINFINITE_SPACE_TR This is an implementation of the Jacobian
%for the semi-infinite space as in REF[1].
%
% NOTATION: REF[1] Section A
%
% - Vector from position a to b: rho_vec_ab = r_b - r_a
% - Distance from a to b: rho_ab = norm(v_rho_ab)
% - Unit vector from position a to b: rho_tilde_ab = rho_vec_ab / rho_ab 
% - sigma = ((gamma + i*omega) / kappa)^ (1/2) for omega > 0
% - sigma_0 = sigma(omega = 0)
% - h_ab = ( ( x_b - x_a ) ^ 2 + ( y_b - y_a ) ^ 2 ) ^ ( 1 / 2 )
% - rho_12_neg = (h_12 ^ 2 + ( z_2 - z_1 ) ^ 2 ) ^ ( 1 / 2 )
% - rho_12_pos = (h_12 ^ 2 + ( z_2 + z_1 ) ^ 2 ) ^ ( 1 / 2 )
% - z_1 = b * ( ( 1 - f_hat ) * mu_s ) ^ (-1) from REF[3]
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


if ~isMellin
    
else
    
end

end

function [J] = jacobian_mellin_moment_semiinfinite_space_TR(z_23, rho_23, rho_12_neg, rho_12_pos)
    

% Compute Jacobian for absorption REF[1] Eq. 47 
J_alpha =  2 * (z_23) / rho_23...
             * (compute_ST(rho_12_neg, rho_23)...
                - compute_ST(rho_12_pos, rho_23));

% Compute Jacobian for diffusion REF[1] Eq. 52
J_nu = - z_12_neg * compute_TT(rho_12_neg, rho_23, omega) ...
       + z_23 * ( rho_tilde_12 * rho_tilde_23 ) * compute_WT(rho_12_neg, rho_23, omega) ...
       - ( - z_12_pos + compute_TT(rho_12_pos, rho23, omega)...
           + z_23 * ( rho_tilde_12 * rho_tilde_23 ) * compute_WT(rho_12_pos, rho_23, omega) );
       

% Return as TODO
J = [J_alpha, J_nu];


end

function [J] = jacobian_temporal_moment_semiinfinite_space_TR()

% Compute Jacobian for absorption REF[1] Eq. 48 
J_alpha =  z_23 / rho_23 * 1 / compute_G(TODO) ...
           * ( compute_ST(rho_12_neg, rho_23) - compute_ST(rho_12_pos, rho_23) ... 
               - compute_t(rho_13) * ( compute_SI(rho_12_neg, rho_23) - compute_ST(rho_12_pos, rho_23) ) );

% Compute Jacobian for diffusion REF[1] Eq. 52
J_nu = 1 / compute_G(TODO)...
       * ( jacobian_mellin_moment_semiinfinite_space_TR(z_23, rho_23, rho_12_neg, rho_12_pos)...
           - compute_t(rho_13) * jacobian_intensity_semiinfinite_space_TR(TODO) );
        

% Return as TODO
J = [J_alpha, J_nu];
end

function S = compute_ST(x, y, kappa, sigma_0)
% Implements REF[1] Eq. 43
    S = 1 / ( 8 * ( 2 * pi ) ^ 3 * kappa ^ 2 * sigma_0)...
        * ( x + sigma_0 * ( x + y ) * y ) / ( x * y ^ 2 )...
        * exp( -sigma_0 * ( x + y ) );
end

function T = compute_TT(x, y, sigma_0, sigma)
% Implements REF[1] Eq. 49.
% TODO where does the hat come from. Eq. 49 defines T but not T_hat
T = 1 / ( 8 * ( 2 * pi ) ^ 3 * kappa ^ 2 )...
    * ( ( x ^ 2 + y ^ 2 ) + sigma_0 * ( x ^ 2 * y + x * y ^ 2 ) ) / ( x * y ) ^ 3 ...
    * exp( - sigma * ( x + y ) );

end

function W = compute_WT(x, y)
% Implements REF[1] Eq. 50.
% TODO where does the hat come from. Eq. 50 defines W but not W_hat
% It might be a typo since earlier in the text (Eq.33) they refer to hat
% symbols as being in the frequency domain, which this isn't.

W = 1 / ( 8 * ( 2 * pi ) ^ 3 * kappa ^ 2 )...
    * ( 3 * ( x ^ 2 + y ^ 2 ) + sigma_0 * ( x * y ^ 2 + 3 * x ^ 2 * y + y ^ 3 ) ...
        + sigma_0 ^ 2 * ( x ^ 2 * y ^ 2 + x * y ^ 3 ) )...
    / ( x ^ 2 * y ^ 4 ); 

end

function t = compute_t(rho, kappa, sigma_0)
% As per REF[1] Eq. 44 and REF[2] Eqs. 4.1 and 4.2 in table 4
% we realise that REF[1] Eq. 44 is also valid for the semi-infinite space
t = 1 / ( 2 * kappa ) * rho ^ 2 / ( 1 + sigma_0 * rho );
end

function G = compute_G(chi_3, r_2)

G = z_23 * ( 1 + sigma * rho_23 ) * exp( - sigma * rho_23)...
    / ( ( 2 * pi ) ^ ( 3 / 2 ) * rho_23 ^ 3 );

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



