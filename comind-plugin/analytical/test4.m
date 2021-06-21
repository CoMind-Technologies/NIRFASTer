clear all, close all

t = 200;            % [ps]
omega = 0;          % [Hz]
order = 1;          % Only moment implemented
c = 0.21;           % [mm/ps]
g = 0.72;           % [dimensionless]
isMellin = true;    % Only one implemented atm  
isRelative = false; % Normalises as in REF[1] Eq. 37
mua = 0.025;        % [mm^-1]
musp = 2.0;         % [mm^-1]
mus = musp / ( 1 - g );

r_1 = [-8.0 0.0 0.0];     % [mm] Source [x, y, z]
r_3 = [ 8.0 0.0 0.0];     % [mm] Detector [x, y, z]

z_depths = (0:50)/4;
x = ((1:100) - 50)/4;

test_label = "CW";
    
J_alpha = zeros(size(z_depths,1),size(x,1));
J_nu = J_alpha;

for ii = 1:length(z_depths)
    disp(z_depths(ii))
    for jj = 1:length(x)
        r_2 = [x(jj), 0, z_depths(ii)]; % [mm] J evaluation [x, y, z]
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
        J_alpha(ii,jj) = J(1);
        J_nu(ii,jj) = J(2);
    end
end

figure('name', 'banana'), subplot(122), hold on
cmin = prctile(J_nu(:), [5]);
cmax = prctile(J_nu(:), [95]);
imagesc(x, z_depths, J_nu, [cmin cmax])
colormap(flipud(hot))
xlabel('x [mm]')
ylabel('z [mm]')
axis('square')
title([ ulabel, '_\nu', '; y[mm] = ', num2str(r_2(2))])
colorbar()

subplot(121), hold on,
cmin = prctile(J_alpha(:), [5]);
cmax = prctile(J_alpha(:), [95]);
imagesc(x, z_depths, J_alpha, [cmin cmax])
colormap(flipud(hot))
xlabel('x [mm]')
ylabel('z [mm]')
axis('square')
title([ ulabel, '_\alpha', '; y[mm] = ', num2str(r_2(2))])
colorbar()

%% Test in NIRFAST

if ~isfile('semiinfinite_slab_fem')
    
    sizevar.width = x(end) - x(1);
    sizevar.height = x(end) - x(1);
    sizevar.depth = z_depths(end) - z_depths(1);
    sizevar.xc = 0;
    sizevar.yc = 0;
    sizevar.zc = (z_depths(end) - z_depths(1))/2;
    sizevar.dist = 1;
    mesh = make_slab(sizevar);

end