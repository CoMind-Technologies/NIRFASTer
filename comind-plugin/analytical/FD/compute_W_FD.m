function W = compute_W_FD(x, y, sigma, kappa)
% Implements REF[1] Eq. Table 2, W^{(\Gamma)}(x,y) for Freq-Domain
W = 1 / ( 4 * ( 2 * pi ) ^ 3 * kappa ) ... % <ok>
    * ( 3 + 3 * sigma * y + sigma ^ 2 * y ^ 2 ) / y ^ 4 ... % <ok>
    * ( 1 + sigma * x ) / x ^ 2 ... % <ok>
    * exp( - sigma * ( x + y ) ); % <ok>
end

