clear all, close all

%% Parameters
t = 200e-9;               % [s]
omega = 0;                % [Hz]
order = 1;                % Only moment implemented
c = 0.21e12;              % [mm/s]
c_vacuum = 0.299792458e12;% [mm/s]
ri = c_vacuum / c;        % refractive index
g = 0.72;                 % [dimensionless]
mua = 0.025;              % [mm^-1]
musp = 2.0;               % [mm^-1]
mus = musp / ( 1 - g );

r_1 = [-8.0 0.0 0.0];     % [mm] Source [x, y, z]
r_3 = [ 8.0 0.0 0.0];     % [mm] Detector [x, y, z]

z_depths = (0:50)/4;
x = ((1:100) - 50)/4;
    
J_alpha = zeros(size(z_depths,1),size(x,1));
J_nu = J_alpha;

test_label = "CW";

%% Test in NIRFAST - crashes if run more than once
isNewMesh = false;   

if isNewMesh
    sizevar.width = x(end) - x(1);
    sizevar.height = z_depths(end) - z_depths(1);
    sizevar.xc = 0;
    sizevar.yc = (z_depths(end) - z_depths(1))/2;
    sizevar.dist = 0.1;
    mesh = make_rectangle(sizevar);

    % keep only one source and one detector
    mesh.source.num = 1;
    mesh.source.coord = [r_1(1), z_depths(end)];
    mesh.source.fwhm = 0;
    mesh.source.fixed = 1;

    mesh.meas.num = 1;
    mesh.meas.fixed = 1;
    mesh.meas.coord = [r_3(1), z_depths(end)];

    mesh.link = [1 1 1];
    mesh.elements = double(mesh.elements);
    mesh.name = "semiinfinite";
    mesh.dimension = 2;

    mesh.type = "stnd";
    mesh.mua = mua * ones(size(mesh.nodes,1),1);
    mesh.mus = musp * ones(size(mesh.nodes,1),1);
    [~, ~, kappa] = compute_optical_parameters(c, mua, mus, g);
    mesh.kappa = kappa / c * ones(size(mesh.nodes,1),1);
    mesh.ri = ri * ones(size(mesh.nodes,1),1);
    mesh.c = c * ones(size(mesh.nodes,1),1); % [mm/s]
    save_mesh(mesh, './rect_mesh')
    mesh = load_mesh('./rect_mesh');
    save('rect_mesh.mat', 'mesh')
end

load('rect_mesh.mat')
plotmesh_fiducials(mesh)

if test_label == "CW"
    [J_FD, ~] = jacobian_FD(mesh, 0, [], [], [], 'all');
else
    error('Not implemented')
end

%% 

plotimage(mesh, 1 / mesh.c(1) * J_FD.complex(1,1:end/2));
title('NIRFAST UNNORM J_{\nu}^{(E)} [mm]');
cl_nu = caxis;

plotimage(mesh, 1 / mesh.c(1) * J_FD.complex(1,(end/2+1):end));
title('NIRFAST UNNORM J_{\alpha}^{(E)} [mm^{-1}]');
cl_alpha = caxis;

plotimage(mesh, J_FD.complete(1,1:end/2));
title('NIRFAST NORM J_{\nu}^{(E)} [mm]');

plotimage(mesh, J_FD.complete(1,(end/2+1):end));
title('NIRFAST NORM J_{\alpha}^{(E)} [mm^{-1}]');


%% Analytical solution

isMellin = false;         % Only one implemented atm  
isRelative = false;        % Normalises as in REF[1] Eq. 37 

for ii = 1:length(z_depths)
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
perc1 = 0;
perc2 = 100;
cmin_nu = prctile(J_nu(:), [perc1]);
cmax_nu = prctile(J_nu(:), [perc2]);
imagesc(x, z_depths, J_nu, cl_nu)
colormap(hot)
xlabel('x [mm]')
ylabel('z [mm]')
axis('square')
if isRelative
    title(['ANALYTIC. NORM ', ulabel, '_\nu', '; y[mm] = ', num2str(r_2(2))])
else
    title(['ANALYTIC. UNNORM ', ulabel, '_\nu', '; y[mm] = ', num2str(r_2(2))])
end
colorbar()

subplot(121), hold on,
cmin_alpha = prctile(J_alpha(:), [perc1]);
cmax_alpha = prctile(J_alpha(:), [perc2]);
imagesc(x, z_depths, J_alpha, cl_alpha)
colormap(hot)
xlabel('x [mm]')
ylabel('z [mm]')
axis('square')
if isRelative
    title(['ANALYTIC. NORM ', ulabel, '_\alpha', '; y[mm] = ', num2str(r_2(2))])
else
     title(['ANALYTIC. UNNORM ', ulabel, '_\alpha', '; y[mm] = ', num2str(r_2(2))])
end
colorbar()

%%

% analytical
ind = find(x == 0);
J_nu_0 = J_nu(:,ind);
J_alpha_0 = J_alpha(:,ind);

figure,
subplot(223), hold on,
plot(z_depths, J_nu_0, 'linewidth', 2)
subplot(221), hold on,
plot(z_depths, J_alpha_0, 'linewidth', 2)

% fem
epsilon = 3e-2;
xinds = find(abs(mesh.nodes(:,1)) <= epsilon);
y_vals = mesh.nodes(xinds,2);
J_nu_1 = 1 / mesh.c(1) * J_FD.complex(1,1:end/2);
J_nu_1 = J_nu_1(1,xinds);
J_nu_1 = interp1(y_vals, J_nu_1,z_depths);

J_alpha_1 = 1 / mesh.c(1) * J_FD.complex(1,end/2:end);
J_alpha_1 = J_alpha_1(1,xinds);
J_alpha_1 = interp1(y_vals, J_alpha_1, z_depths);

subplot(223), hold on,
plot(z_depths, J_nu_1, 'linewidth', 2)
legend('Analytical', 'NIRFAST')
ylabel('J_{\nu}')
xlabel('z [mm]')

subplot(221), hold on,
plot(z_depths, J_alpha_1, 'linewidth', 2)
legend('Analytical', 'NIRFAST')
ylabel('J_{\alpha}')
xlabel('z [mm]')
title('NON SCALED')

% norm
J_nu_0 = -J_nu_0;
mn = min(J_nu_0);
mx = max(J_nu_0);
J_nu_0_norm = (J_nu_0 - mn) / (mx - mn);
mn = min(J_alpha_0);
mx = max(J_alpha_0);
J_alpha_0_norm = (J_alpha_0 - mn) / (mx - mn);

mn = min(J_nu_1);
mx = max(J_nu_1);
J_nu_1_norm = (J_nu_1 - mn) / (mx - mn);

mn = min(J_alpha_1);
mx = max(J_alpha_1);
J_alpha_1_norm = (J_alpha_1 - mn) / (mx - mn);

subplot(224), hold on,
plot(z_depths, J_nu_0_norm, 'linewidth', 2)
subplot(222), hold on,
plot(z_depths, J_alpha_0_norm, 'linewidth', 2)
title('SCALED')

subplot(224), hold on,
plot(z_depths, J_nu_1_norm, 'linewidth', 2)
xlabel('z [mm]')
ylabel('J_{\nu}')
legend('Analytical', 'NIRFAST')

subplot(222), hold on,
plot(z_depths, J_alpha_1_norm, 'linewidth', 2)
xlabel('z [mm]')
ylabel('J_{\alpha}')
legend('Analytical', 'NIRFAST')


%% norm

% analytical
ind = find(x == 0);
J_nu_0 = J_nu(:,ind);
J_alpha_0 = J_alpha(:,ind);

figure,
subplot(223), hold on,
plot(z_depths, J_nu_0, 'linewidth', 2)
subplot(221), hold on,
plot(z_depths, J_alpha_0, 'linewidth', 2)

% fem
epsilon = 3e-2;
xinds = find(abs(mesh.nodes(:,1)) <= epsilon);
y_vals = mesh.nodes(xinds,2);
J_nu_1 = J_FD.complete(1,1:end/2);
J_nu_1 = J_nu_1(1,xinds);
J_nu_1 = interp1(y_vals, J_nu_1,z_depths);

J_alpha_1 = J_FD.complete(1,end/2:end);
J_alpha_1 = J_alpha_1(1,xinds);
J_alpha_1 = interp1(y_vals, J_alpha_1, z_depths);

subplot(223), hold on,
plot(z_depths, J_nu_1, 'linewidth', 2)
legend('Analytical', 'NIRFAST')
ylabel('J_{\nu}')
xlabel('z [mm]')

subplot(221), hold on,
plot(z_depths, J_alpha_1, 'linewidth', 2)
legend('Analytical', 'NIRFAST')
ylabel('J_{\alpha}')
xlabel('z [mm]')
title('NON SCALED')

% norm
J_nu_0 = -J_nu_0;
mn = min(J_nu_0);
mx = max(J_nu_0);
J_nu_0_norm = (J_nu_0 - mn) / (mx - mn);
mn = min(J_alpha_0);
mx = max(J_alpha_0);
J_alpha_0_norm = (J_alpha_0 - mn) / (mx - mn);

mn = min(J_nu_1);
mx = max(J_nu_1);
J_nu_1_norm = (J_nu_1 - mn) / (mx - mn);

mn = min(J_alpha_1);
mx = max(J_alpha_1);
J_alpha_1_norm = (J_alpha_1 - mn) / (mx - mn);

subplot(224), hold on,
plot(z_depths, J_nu_0_norm, 'linewidth', 2)
subplot(222), hold on,
plot(z_depths, J_alpha_0_norm, 'linewidth', 2)
title('SCALED')

subplot(224), hold on,
plot(z_depths, J_nu_1_norm, 'linewidth', 2)
xlabel('z [mm]')
ylabel('J_{\nu}')
legend('Analytical', 'NIRFAST')

subplot(222), hold on,
plot(z_depths, J_alpha_1_norm, 'linewidth', 2)
xlabel('z [mm]')
ylabel('J_{\alpha}')
legend('Analytical', 'NIRFAST')



