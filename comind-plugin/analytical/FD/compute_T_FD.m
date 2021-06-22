function T = compute_T_FD(x, y, sigma, kappa) % <ok>
% Implements REF[1] Eq. Table 2, T^{(\Gamma)}(x,y) for Freq-Domain
T = 1 / ( 4 * ( 2 * pi ) ^ 3 * kappa ) ... % <ok>
    * ( 1 + sigma * y ) / y ^ 3 ... % <ok>
    * ( 1 + sigma * x ) / x ^ 3 ... % <ok>
    * exp( - sigma * ( x + y ) );   % <ok>
end
