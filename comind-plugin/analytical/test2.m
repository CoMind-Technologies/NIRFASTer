clear all, close all

t = 200;            % [ps]
omega = 0;          % [Hz]
order = 1;          % Only moment implemented
c = 0.21;           % [mm/ps]
g = 0.72;           % [dimensionless]
isMellin = true;   % Only one implemented atm
isRelative = false;
mua = 0.025;        % [mm^-1]
musp = 2.0;         % [mm^-1]
mus = musp / ( 1 - g );
z_depth = 3;

test_label = "T";


    
r_1 = [-8.0 0.0 0.0];     % [mm] source [x, y, z]
r_2 = [ 0.0 0.0 z_depth]; % [mm] J evaluation [x, y, z]
r_3 = [ 8.0 0.0 0.0];     % [mm] detector [x, y, z]


x = ((1:100) - 50)/4;
J_alpha = zeros(size(x,1),size(x,1));
J_nu = J_alpha;

for i = 1:length(x)
    for j = 1:length(x)
        r_2 = [x(j), x(i), r_2(3)];
        if test_label == "CW"
            J = jacobian_integrated_intensity_semiinfinite_space_CW(r_1, r_2, r_3, mua, c, mus, g, isRelative);
            ulabel = 'J^{(E)}';
        elseif test_label == "flux"
            J = jacobian_flux_time_semiinfinite_space_TD(r_1, r_2, r_3, t, mua, c, mus, g); %<ok>
            ulabel = ['J^{(\Gamma)}(t[ps] = ', num2str(t), ')'];
        elseif test_label == "T"
            J = jacobian_moment_semiinfinite_space_FD(r_1, r_2, r_3, mua, c, mus, g, order, isMellin);
            ulabel = 'J^{(T)}';
        else
            error("Not implemented")
        end
        J_alpha(i,j) = J(1);
        J_nu(i,j) = J(2);
    end
end

figure('name', test_label), subplot(1,2,2), hold on,
[X, Y] = meshgrid(x);
plot3(r_1(1), r_1(2), r_1(3), '*', 'linewidth', 3)
plot3(r_3(1), r_3(2), r_3(3), 'r^', 'linewidth', 3)
mesh(X, Y, J_nu, 'linewidth', 2)
zlabel([ulabel,'_\nu'])
xlabel('x [mm]')
ylabel('y [mm]')
axis('square')
title(['z = ', num2str(z_depth)])
view(20,26)
grid()

subplot(1,2,1), hold on,
plot3(r_1(1), r_1(2), r_1(3), '*', 'linewidth', 3)
plot3(r_3(1), r_3(2), r_3(3), 'r^', 'linewidth', 3)
mesh(X, Y, J_alpha, 'linewidth', 2)
zlabel([ulabel,'_\alpha'])
xlabel('x [mm]')
ylabel('y [mm]')
axis('square')
title(['z = ', num2str(z_depth)])
view(20,26)
grid()

