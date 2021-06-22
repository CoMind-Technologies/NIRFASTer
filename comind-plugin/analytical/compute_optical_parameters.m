function [z_1, gamma, kappa] = compute_optical_parameters(c, mua, mus, g)
    musp = ( 1 - g ) * mus;
    z_1 = 1 / musp;
    gamma = mua * c;
    kappa = c / ( 3 * ( mua + musp ) );
end

