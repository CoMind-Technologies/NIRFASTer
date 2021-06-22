function [rho_23, rho_12_neg, rho_12_pos, rho_tilde_12_neg, ...
    rho_tilde_12_pos, rho_tilde_23, z_23, z_12_neg, z_12_pos] = compute_positions(r_1, r_2, r_3, z_1)
    r_1_neg = r_1 + [0 0 z_1]; % By convention in the z+ side
    r_1_pos = r_1 - [0 0 z_1]; % By convention in the z- side (mirror source)
    rho_23 = norm(r_3 - r_2);
    rho_12_neg = compute_rho_ab(r_1_neg, r_2);
    rho_12_pos = compute_rho_ab(r_1_pos, r_2);

    rho_tilde_12_neg = compute_rho_ab_tilde(r_1_neg, r_2);
    rho_tilde_12_pos = compute_rho_ab_tilde(r_1_pos, r_2);
    rho_tilde_23 = compute_rho_ab_tilde(r_2, r_3);

    % For a detector s(r_3) = [0, 0, -1], this are the z-dim of the equivalent
    % r_ab = r_b - r_a vectors
    z_23 = (r_3(3) - r_2(3)); % <ok>
    z_12_neg = (r_2(3) - r_1_neg(3)); % <ok>
    z_12_pos = (r_2(3) - r_1_pos(3)); % <ok>
end

function rho_tilde = compute_rho_ab_tilde(r_a, r_b)
    rho_tilde = ( r_b - r_a ) / norm( r_b - r_a );
end

function h = compute_h_ab(r_a, r_b)
    % Implements REF[1] Table 1 row 2 and 3
    % h = ( ( x_b - x_a ) ^ 2 + ( y_b - y_a ) ^ 2 ) ^ ( 1 / 2 );
    h = ( ( r_b(1) - r_a(1) ) ^ 2 + ( r_b(2) - r_a(2) ) ^ 2 ) ^ ( 1 / 2 );
end

function rho = compute_rho_ab(r_a, r_b)
    % Implements REF[1] Table 1 row 2 and 3
    % rho_ab_neg = (h_ab ^ 2 + ( z_b - z_a ) ^ 2 ) ^ ( 1 / 2 )
    h_ab = compute_h_ab(r_a, r_b);
    rho = ( h_ab ^ 2 + ( r_b(3) - r_a(3) ) ^ 2 ) ^ ( 1 / 2 );
end

