clear all, close all,

t = 200;            % [ps]
c = 0.21;           % [mm/ps]
g = 0.72;           % [dimensionless]
isMellin = true;    % Only one implemented atm        
mua = 0.025;        % [mm^-1]
musp = 2.0;         % [mm^-1]
mus = musp / ( 1 - g );

for k = 1:1
    
    r_1 = [-8.0 0.0 0.0]; % [mm] source [x, y, z]
    r_2 = [ 0.0 0.0 0.9]; % [mm] J evaluation [x, y, z]
    r_3 = [ 8.0 0.0 0.0]; % [mm] detector [x, y, z]
    
    
    x = ((1:100) - 50)/4;
    J_alpha = zeros(size(x,1),size(x,1));
    J_nu = J_alpha;
    
    for i = 1:length(x)
        for j = 1:length(x)
            r_2 = [x(i), x(j), r_2(3)];
            J = jacobian_flux_time_semiinfinite_space_TR(r_1, r_2, r_3, t, mua, c, mus, g);
            J_alpha(i,j) = J(1);
            J_nu(i,j) = J(2);
        end
    end

    figure, hold on,
    [X, Y] = meshgrid(x);
    plot3(r_1(1), r_1(2), r_1(3), '*', 'linewidth', 3)
    plot3(r_3(1), r_3(2), r_3(3), 'r^', 'linewidth', 3)
    mesh(X', Y', J_nu, 'linewidth', 2)
    zlabel('J^{(\Gamma)}_\nu')
    xlabel('x [mm]')
    ylabel('y [mm]')
    view(20,26)
    
    figure, hold on,
    plot3(r_1(1), r_1(2), r_1(3), '*', 'linewidth', 3)
    plot3(r_3(1), r_3(2), r_3(3), 'r^', 'linewidth', 3)
    mesh(X', Y', J_alpha, 'linewidth', 2)
    zlabel('J^{(\Gamma)}_\alpha')
    xlabel('x [mm]')
    ylabel('y [mm]')
    view(20,26)

end