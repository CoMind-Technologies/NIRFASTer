function S = compute_S_TD(x, y, t, kappa, gamma) %<ok>
    S = 1 / ( 2 * ( 2 * pi ) ^ 3 )... %<ok>
        * ( 1 / y ^ 2 + ( x + y ) / ( 2 * kappa * t ) * ( 1 / x + 1 / y ) )... %<ok>
        * compute_f_TD(x + y, t, kappa, gamma); %<ok>
end
