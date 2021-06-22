function [node_errors] = fem_analytical_jacobians_error(homog_mesh,...
                                       jac_analytical_handle,...
                                       jac_fem_handle,...
                                       jac_analytical_cfg,...
                                       jac_fem_cfg)
                                   
% RETURNS
% NODE_ERRORS: [n_node, 1] of reconstruction errors at the location of each 
%   mesh node 
%
% ARGS
% HOMOG_MESH: a NIRFAST mesh which is homogenous in optical properties and
%   which has just one source and one detector
% JAC_ANALYTICAL_HANDLE: a function handle for the function that computes 
%   the analytical Jacobian 
% JAC_FEM_HANDLE: a handle for the FEM Jacobian function we want to test
% JAC_ANALYTICAL_CFG: a struct of fields and values to set the optional 
%   argument values in the analytical jacobian function (see the
%   particular function for details of arguments)
% JAC_FEM_CFG: a struct of fields and values to set the optical arguement
%   values in the fem jacobna function (see the particular function for
%   details of arguments)


% TODO:
% * In a wrapper script, we need to generate the homogenous mesh in order to
%  pass it in here. This mesh contains info on the sd positions
%  (separations) and therefore needs to change for a range of resutls (e.g.
%  curve of approx error against sd separation). So first build mesh, then
%  simply modify the sd positions field in place before passing to this
%  function
% * (Pablo) write function which returns Jacobian based on analytical
%  solution, mimicking the function signature of
%  NIRFASTer.analytical.semi_infinite_TR
% * fem_TR jacobian function and analytical_TR jacobian function should be
%  supplied with the same T and DT

% checks
if (length(unique(homog_mesh.mua)) ~= 1) || (length(unique(homog_mesh.mua)) ~= 1) || (length(unique(homog_mesh.ri)) ~= 1)
    error('Mesh supplied must be homogenous in optical properties');
else
    mu_a = homog_mesh.mua(1);
    mu_sp = homog_mesh.mus(1);
    ri = homog_mesh.ri(1);
end
if (numel(homog_mesh.source.num) ~= 1)|| (numel(homog_mesh.meas.num) ~= 1)
    error('Mesh must have just one source and one detector');
end

% derive some constants
c0 = 299792458000; %/mm*s-1 - speed of light in air/vaccum
cm = c0/ri; %/mm*s-1 - speed of light in medium

% unpack configs
allowed_args_fem = ["order", "central", "second_mesh_basis", "solver", "options", "all"];
fem_cfg_names = fieldnames(jac_fem_cfg);
fem_parsed_arg_cell = {};
for ni = 1:length(allowed_args_fem)
    nn = allowed_args_fem(ni);
    if any(strcmp(nn, fem_cfg_names))
        if strcmp(ann, 'central')
            fem_parsed_arg_list{ni} = 'central';
        else
            fem_parsed_arg_lsit{ni} = getfield(jac_fem_cfg, nn);
        end
    else
        fem_parsed_arg_list{ni} = [];
    end
end

% We must have an anisoptropy supplied in addition to the NIRFAST mesh
% since the analytical jacobian computation requires mu_s; not mu_s'
if any(strcmp(fem_cfg_names, 'g'))
    g = jac_fem_cfg.g;
elseif any(strcmp(fem_cfg_names, 'anisotropy'))
    g = jac_fem_cfg.anisotropy;
else
    error('You must supply an anistropy in the fem config');
end
mu_s = mu_sp / (1-g);

allowed_args_an = ["omega", "isMellin"];
an_cfg_names = fieldnames(jac_analytical_cfg);
if any(strcmp(an_cfg_names, 'omega'))
    omega = jac_analytical_cfg.omega;
else
    omega = 0;  % default
end
if any(strcmp(an_cfg_names, 'isMellin'))
    isMellin = jac_analytical_cfg.isMellin;
else
    isMellin = false;  % default
end

% an_parsed_arg_cell = {};
% for ni = 1:length(allowed_args_an)
%     nn = allowed_args_an(ni);
%     if any(strcmp(nn, an_cfg_names))
%         % TODO
%         an_parsed_arg_cell{ni} = getfield(jac_analytical_cfg, nn);
%     else
%         an_parsed_arg_cell{ni} = [];
%     end
% end

% call the fem jacobian function
J_fem = jac_fem_handle(homog_mesh, fem_parsed_arg_cell{:});
% [kappa_1,..., kappa_N, mu_a_1,..., mu_a_N] % (1)
J_fem = reshape(J_fem, numel(J_fem)/2, 2)';

source_coords = homog_mesh.source.coord;
det_coords = homog_mesh.meas.coord;
sd_distance = norm(source_coords - det_coords);

% evaluate the analytical Jacobian at the sd distance given by the mesh's 
% source and detector locations% and at the locations of only the mesh 
% nodes (so we can compare like for like)
node_errors = zeros(size(homog_mesh.nodes, 1), 1);
for ni = 1:size(homog_mesh.nodes, 1)
    node_coords = homog_mesh.nodes(ni,:);
    analytical_node = jac_analytical_handle(source_coords, node_coords,...
                                            det_coords, omega, mu_a, c,...
                                            mu_s, g, isMellin);
    node_errors(ni) = fliplr(analytical_node) - J_fem(:, ni)';
end
    






