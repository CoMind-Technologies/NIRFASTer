function W = compute_W_TD(x, y, t, kappa, gamma) %<ok>
% Implements REF[1] Eq. Table 2, W^{(\Gamma)}(x,y) for Time-Domain
W = (...  %<ok>
        3 * x / ( y ^ 4 * t ) ...  %<ok>
        + ( x + y ) ^ 2 ...  %<ok>
          * ( 3 * x ^ 2 - 2 * x * y + y ^ 2 ) ...  %<ok>
          / ( x ^ 2 * y ^ 3 ) ...  %<ok>
          * ( 1 / ( 2 * kappa * t ^ 2 ) ) ...  %<ok>
        + ( x + y ) ^ 4 / ( x * y ^ 2) * 1 / ( 4 * kappa ^ 2 * t ^ 3) ...  %<ok> 
    )...  %<ok>
    * 1 / ( 4 * kappa * ( 2 * pi ) ^ 3 )...  %<ok>
    * compute_f_TD(x + y, t, kappa, gamma);   %<ok>
end  %<ok>
