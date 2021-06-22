function T = compute_T_TD(x, y, t, kappa, gamma)  %<ok>
% Implements REF[1] Eq. Table 2, T^{(\Gamma)}(x,y) for Time-Domain
T = (...  %<ok>
        (1 / x  + 1 / y ) ^ 3 ...  %<ok>
          * ( x * y / ( 2 * kappa * t ^ 2 ) ) ...  %<ok>
        + ( 1 / x ^ 3  + 1 / y ^ 3 )...  %<ok>
          * 1 / t ...  %<ok>
    )...  %<ok>
    * 1 / ( 4 * kappa * ( 2 * pi ) ^ 3 )...  %<ok>
    * compute_f_TD(x + y, t, kappa, gamma);  %<ok>
end
