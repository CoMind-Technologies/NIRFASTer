function f = compute_f_TD(x, t, kappa, gamma) %<ok>
% Implements REF[1] Eq. Table 2, T^{(\Gamma)}(x,y) for Time-Domain
f = ( 2 * kappa * t ) ^ ( - 3 / 2 )... %<ok>
    * exp( - gamma * t - x ^ 2 / ( 4 * kappa * t ) ); %<ok>
end

