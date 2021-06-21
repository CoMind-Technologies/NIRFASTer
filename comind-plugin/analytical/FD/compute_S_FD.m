function S = compute_S_FD(x, y, sigma, kappa) % <ok>
% Implements REF[1] Eq. Table 2, S^{(\Gamma)}(x,y) for Freq-Domain
S = 1 / ( 4 * ( 2 * pi ) ^ 3 * kappa )... % <ok>
    * ( 1 + sigma * y ) / ( x * y ^ 2 ) ... % <ok>
    * exp( - sigma * ( x + y ) ); % <ok>
end

