function J = jacobian_integrated_intensity_semiinfinite_space_CW(r_1, r_2, r_3, mua, c, mus, g, isRelative)
% This is an implementation of the Jacobian
%for the semi-infinite space as in REF[1].
%
% REFS:
%   [1] Arridge, S. R. (1995). Photon-measurement density functions. 
%       Part I: Analytical forms. Applied Optics, 34(31), 7395-7409.

% REF[1] Section B. " The jacobian for the integrated intenstiy J^{(E)} is trivially obtained 
% as the DC component of the frequency domain jacobian for the flux  J(omega=0)^{(\Gamma)}"

omega = 0; % <ok>
J = jacobian_flux_semiinfinite_space_FD(r_1, r_2, r_3, mua, c, mus, g, omega); % <ok>



if isRelative
    [~, gamma, kappa] = compute_optical_parameters(c, mua, mus, g);
    sigma_0 = compute_sigma(gamma, 0, kappa);
    G = compute_G_FD(r_3, r_2, sigma_0);
    J = J / G;
end

end

