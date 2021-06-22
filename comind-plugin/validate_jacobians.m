
szv.width=30;
szv.height=30;
szv.depth=60;
szv.xc=0;
szv.yc=0;
szv.zc=0;
szv.dist=0.7;
mesh_semi = make_slab(szv);
mesh_semi.elements = double(mesh_semi.elements);

mu_a = 0.025;
mu_sp = 2.0;
c0 = 299792458000;
c = 0.21;
source_coords = [-8, 0, -30];
det_coords = [8, 0, -30];

mesh_semi.mua = mu_a .* ones(size(mesh_semi.nodes, 1), 1);
mesh_semi.mus = mu_sp .* ones(size(mesh_semi.nodes, 1), 1);
mesh_semi.kappa = ones(size(mesh_semi.mua)) ./ (3 * (mesh_semi.mua + mesh_semi.mus));
mesh_semi.ri = (c0/c) .* ones(size(mesh_semi.nodes, 1), 1);
mesh_semi.c = c .* ones(size(mesh_semi.nodes, 1), 1);

%%%
mesh_semi.ksi = 0.2129 .* ones(size(mesh_semi.nodes, 1), 1);
%%%

mesh_semi.source.fixed = 0;
mesh_semi.source.num = 1;
mesh_semi.source.coord = source_coords;
mesh_semi.source.fwhm = 0;
mesh_semi.source.distributed = 0;

mesh_semi.meas.fixed = 1;
mesh_semi.meas.num = 1;
mesh_semi.meas.coord = det_coords;
mesh_semi.link = [1, 1, 1];

% NB currently need to save and then (re)load mesh for int_func to
% magically appear!

[J_fem, ~] = jacobian_stnd_FD(mesh_semi, 0, [], [], [], 'all');  % CW
% [J_fem, ~] = jacobian_stnd_TR(mesh_semi, 1, [], [], [], [], 'all');
J_fem_ = J_fem.first.complete;
J_kappa = J_fem_(1:numel(J_fem_)/2);
J_mua = J_fem_(numel(J_fem_)/2+1:end);
J_mus = (c.*ones(size(J_kappa))./3.*J_kappa) - J_mua;

[x,y,z]=meshgrid(min(mesh_semi.nodes(:,1)):1:max(mesh_semi.nodes(:,1)),15,min(mesh_semi.nodes(:,3)):1:max(mesh_semi.nodes(:,3)));
F = TriScatteredInterp(mesh_semi.nodes(:,1),mesh_semi.nodes(:,2),mesh_semi.nodes(:,3),J_mua');
val = F(x,y,z);
figure; surf(squeeze(x),squeeze(y),squeeze(z),squeeze(val));
view(0,0); shading interp; colorbar; axis tight

[x,y,z]=meshgrid(0,0,min(mesh_semi.nodes(:,3)):1:max(mesh_semi.nodes(:,3)));
F = TriScatteredInterp(mesh_semi.nodes(:,1),mesh_semi.nodes(:,2),mesh_semi.nodes(:,3),J_mus');
val = F(x,y,z);
figure; plot(squeeze(val))

figure, hold on,
% [X, Y] = meshgrid(x);
[x,y, z]=meshgrid(min(mesh_semi.nodes(:,1)):1:max(mesh_semi.nodes(:,1)),min(mesh_semi.nodes(:,2)):1:max(mesh_semi.nodes(:,2)), 1);
% plot3(r_1(1), r_1(2), r_1(3), '*', 'linewidth', 3)
% plot3(r_3(1), r_3(2), r_3(3), 'r^', 'linewidth', 3)
mesh(X', Y', J_nu, 'linewidth', 2)
zlabel('J^{(\Gamma)}_\nu')
xlabel('x [mm]')
ylabel('y [mm]')
view(20,26)


