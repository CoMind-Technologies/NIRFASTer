function t_Gamma = compute_t_Gamma_FD(x, sigma, kappa)
% As per REF[1] Eq. 44 and REF[2] Eqs. 4.1 and 4.2 in table 4
% we realise that REF[1] Eq. 44 is also valid for the semi-infinite space
t_Gamma =  1 / ( 2 * kappa ) * x ^ 2 / ( 1 + sigma * x ); % <ok>
end
