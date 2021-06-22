function G = compute_G_FD(r_a, r_b, sigma)
% REFS:
%   [1] Arridge, S. R. (1995). Photon-measurement density functions. 
%       Part I: Analytical forms. Applied Optics, 34(31), 7395-7409.

rho_ab = r_b - r_a;

% REF[1] Eq. 26
G = rho_ab(3) * ( 1 + sigma * norm(rho_ab) ) * exp( - sigma * norm(rho_ab) )...
    / ( ( 2 * pi ) ^ ( 3 / 2 ) * norm(rho_ab) ^ 3 );

end

