function sigma = compute_sigma(gamma, omega, kappa)
    % Implements REF[1] Table 1 row 1
    sigma = ((gamma + 1i*omega) / kappa) ^ (1/2);
end