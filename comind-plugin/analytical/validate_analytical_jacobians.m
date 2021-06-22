function [error] = validate_analytical_jacobians(jac_analytical_handle,...
                                                 mu_a_0,...
                                                 mu_s_0,...
                                                 delta_mu_a,...
                                                 delta_mu_s,...
                                                 ri,...
                                                 T,...
                                                 dt)

% Validate the analytical derivation of jacobians (for mu_a and kappa) by
% comparing against delta_phi / delta_mu numerically, where delta_mu is a
% supplied delta in [kappa, mu_a] from the initial values specified in
% mesh, and delta_phi is phi(mu_0 + delta_mu) - phi(mu_0) as computed by
% NIRFASTer's analytical semi-infinite solution

% TODO:
% * see warning in semi_infinite_tr - what's the significance of this?
% * so we want dphi or dgamma where gamma is boundary measurement?
% * use semi_infinite_TD or semi_infinite_FD?

phi_0_all_dist = [];
phi_1_all_distm= [];
for sd_distance = 1:25 % CHECK: mm? is this a sensible range?
    phi_0 = semi_infinite_TR(mu_a_0, mu_s_0, ri, sd_distance, T, dt);
    phi_0_all_dist = [phi_0_all_dist, phi_0];
    phi_1 = semi_infinite_TR(mu_a_0+delta_mu_a, mu_s_0+delta_mu_s, ri, sd_distance, T, dt);
    phi_1_all_dist = [phi_1_all_dist, phi_1];
    
end

% TODO: convert to non-TR or keep as TR??? Assume former for now...


error = repmat(phi_1 - phi_0, 1, 2) ./ [delta_mu_a, delta_mu_s];

