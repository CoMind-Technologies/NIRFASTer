% Try to recreate the plots of FEM Jacobian along a 45 degree radial line of
% a circle in "Photon-measurement density functions. Part 2:
% Finite-element-method calculations", Arridge and Schweiger 1995

% NB we use the built-in Nirfast 2D circle which has different number of
% nodes and elements to the one in the paper, but the plots are normalised
% and so should be closely comparable. The paper mesh also seems to have a
% double mesh density at the boundary, while nf mesh does not. This could
% be more of a problem for comparison.

mesh_c = load_mesh('circle2000_86_stnd');

mu_sp = 2;
mu_a_list = [0.0025, 0.005, 0.01, 0.02, 0.04, 0.08];


chosen_source = 13;  % nearest to 12 o'clock
chosen_det = 1; % nearest to 3 o'clock
chosen_link = find(sum(mesh_c.link(:, 1:2)==[chosen_source, chosen_det], 2) == 2);

mesh_c.source.num = mesh_c.source.num(1,:);  % stupid hack
mesh_c.source.coord = mesh_c.source.coord(chosen_source,:);
mesh_c.source.fwhm = mesh_c.source.fwhm(chosen_source,:);
mesh_c.source.int_func = mesh_c.source.int_func(chosen_source,:);

mesh_c.meas.num = mesh_c.meas.num(1,:);  % stupid hack
mesh_c.meas.coord = mesh_c.meas.coord(chosen_det,:);
mesh_c.meas.int_func = mesh_c.meas.int_func(chosen_det,:);

mesh_c.link = [1 1 1];  % stupid hack

% COMMENTED OUT BELOW is what I think should be the way to simulate only
% the selected pairs, but it crashes `build_jacobian_FD` so for now, sim
% all of the pairs and then just recover the pairwise sub-J of interest

% on_bool = (sum(mesh_c.link == [chosen_source, chosen_det, 1], 2)) == 3;
% mesh_c.link(:, 3) = on_bool;


mesh_c.mus = ones(size(mesh_c.mus)) .* mu_sp;

% I think in the paper, J is normalised by ri, since no ri value for the
% mesh is given. So hopefully if we set it to 1 then we don't need to
% adjust the Nirfast J
mesh_c.ri = ones(size(mesh_c.ri)) .* (0.21e12/0.279e12);  % check this
mesh_c.c = ones(size(mesh_c.c)) .* 0.21e12;

% NB I have no idea what to do about ksi!!!

generate_plots(mesh_c, mu_a_list, mu_sp, chosen_source, chosen_det, chosen_link, 'zeroth', false, true)



function generate_plots(mesh_c, mu_a_list, mu_s, chosen_source, chosen_det, chosen_link, jac_type, normalised, normalised1)

if normalised
    plot_title_arg1 = '(normalised)';
else
    plot_title_arg1 = '(unnormalised)'; 
end
if normalised1
    plot_title_arg2 = '(max=1-normalised)';
else
    plot_title_arg2 = '';
end

f1 = figure;

f2 = figure;


for mi = 1:length(mu_a_list)
    mu_a = mu_a_list(mi);
    mesh_c.mua = ones(size(mesh_c.mua)) .* mu_a;
    mesh_c.mus = ones(size(mesh_c.mus)) .* mu_s;
    mesh_c.kappa = ones(size(mesh_c.mua)) ./ (3.*(mesh_c.mua + mesh_c.mus));

    % which Jacobian type
    if any(strcmp(jac_type, ["zeroth", "fd", "amplitude"]))
        [J_fem, ~] = jacobian_stnd_FD(mesh_c, 0, [], [], [], 'all'); % NB this is dlogI/dmu
    elseif any(strcmp(jac_type, ["first", "mean"]))
        [J_fem, ~] = jacobian_stnd_TR(mesh_c, 1, [], [], [], [], 'all');
        J_fem = J_fem.first;
    elseif any(strcmp(jac_type, ["second", "variance"]))
        [J_fem, ~] = jacobian_stnd_TR(mesh_c, 2, [], [], [], [], 'all');
        J_fem = J_fem.second;
    end
    if normalised
    %     plot_title_arg = '(normalised)';
        J_fem = J_fem.complete;
    else
    %     plot_title_arg = '(unnormalised)';
        J_fem = J_fem.complex;
    end

%     J_fem= J_fem(chosen_link, :);
    % separate mu_a and kappa parts of J
    J_kappa = J_fem(1:numel(J_fem)/2);
    J_mua = J_fem(numel(J_fem)/2+1:end);
  

    % plotimage(mesh_c, J_mua);
    % title('Single SD Jacobian: $\mu_a$', jac_type, plot_title_arg, 'Interpreter','latex');
    % 
    % plotimage(mesh_c, J_kappa);
    % title('Single SD Jacobian: $\mu_{\kappa}$', jac_type, plot_title_arg, 'Interpreter','latex');

    % interpolate onto S-D bisecting radial line
    % first get gradient of radial line from gradient of source-det chord
%     m_chord = (mesh_c.source.coord(chosen_source, 2)-mesh_c.meas.coord(chosen_det, 2))/(mesh_c.source.coord(chosen_source, 1)-mesh_c.meas.coord(chosen_det, 1));
    m_chord = (mesh_c.source.coord(1, 2)-mesh_c.meas.coord(1, 2))/(mesh_c.source.coord(1, 1)-mesh_c.meas.coord(1, 1));
    m_rad = -1.0/m_chord;
    % do the interpolation via a meshgrid
    X = meshgrid(0:0.1:max(mesh_c.nodes(:,1)));
    %     [X, Y] = meshgrid(min(mesh_c.nodes(:,1)):0.5:max(mesh_c.nodes(:,1)), min(mesh_c.nodes(:,2)):0.5:max(mesh_c.nodes(:,2)));
    Y = m_rad .* X;
    F_mua = scatteredInterpolant(mesh_c.nodes(:,1), mesh_c.nodes(:,2), J_mua');
    F_kappa = scatteredInterpolant(mesh_c.nodes(:,1), mesh_c.nodes(:,2), J_kappa');
    J_rad_mua = F_mua(X, Y);
    J_rad_kappa = F_kappa(X, Y);
    % all rows are repeated due to meshgrid, so just take first
    J_rad_mua = J_rad_mua(1,:);
    J_rad_kappa = J_rad_kappa(1,:);
    % restrict grid points to lie strictly within the circle volume
    x = X(1,:);
    y = Y(1,:);
    circle_radius = min(max(mesh_c.nodes(:,1)), max(mesh_c.nodes(:,2))); % slighlty wonky
    for i = 1:numel(x)
        xi = x(i);
        yi = y(i);
        if norm([xi, yi]) >= circle_radius
            break
        end
    end
    x = x(1:i);
    y = y(1:i);
    J_rad_mua = J_rad_mua(1:i);
    J_rad_kappa = J_rad_kappa(1:i);

    % radial distance from centre

    xy = [x; y];
    r = vecnorm(xy);
    
    % possibly normalise so that max (min) value along radius is 1 (-1)
    if normalised1
        maxkappa = max(J_rad_kappa);
        if maxkappa > 0
            J_rad_kappa = J_rad_kappa ./ maxkappa;
        else
            J_rad_kappa = J_rad_kappa ./ min(J_rad_kappa);
        end
        maxmua = max(J_rad_mua);
        if maxmua > 0
            J_rad_mua = J_rad_mua ./ maxmua;
        else
            J_rad_mua = J_rad_mua ./ min(J_rad_mua);
        end
    end

    % plot
    figure(f1); plot(r, J_rad_mua, 'DisplayName', ['mu_a=', num2str(mu_a)]);
    hold on;
    % title('Single SD Jacobian along SD bisecting radius: $\mu_a$', jac_type, plot_title_arg, 'Interpreter','latex');
    % xlabel('distance from centre [mm]')
    figure(f2); plot(r, J_rad_kappa, 'DisplayName', ['mu_a=', num2str(mu_a)]);
    hold on;
    % title('Single SD Jacobian along SD bisecting radius: $\mu_{\kappa}$', jac_type, plot_title_arg, 'Interpreter','latex');
    % xlabel('distance from centre [mm]')
end

figure(f1); legend;
title('Single SD Jacobian along SD bisecting radius: $\mu_a$', [jac_type, ' ', plot_title_arg1, ' moment', plot_title_arg2], 'Interpreter','latex');
xlabel('distance from centre [mm]')
figure(f2); legend;
title('Single SD Jacobian along SD bisecting radius: $\mu_{\kappa}$', [jac_type, ' ', plot_title_arg1, ' moment', plot_title_arg2], 'Interpreter','latex');
xlabel('distance from centre [mm]')

end

    
    
       
    
    
