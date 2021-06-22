clear all, close all

omega = 0;       % [Hz] We want to evaluate moments in the stationary condition
                 % i.e. when the full TPSF has been captured
                 % i.e. when Lambda(t) when t --> inf
                 % which in Fourier space corresponds to Lambda(omega = 0)
order = 1;
t = 200;         % [ps]
c = 0.21;        % [mm/ps]
g = 0.72;        % [dimensionless]
isMellin = false; % Only one implemented atm  
isRelative = false
mua = 0.025;    % [mm^-1]
musp = 2.0;     % [mm^-1]
mus = musp / ( 1 - g );
num_z = 1000;

figure,

r_1 = [-8.0 0.0 0.0]; % [mm] source [x, y, z]
r_2 = [ 0.0 0.0 0.0]; % [mm] J evaluation [x, y, z]
r_3 = [ 8.0 0.0 0.0]; % [mm] detector [x, y, z]

JT = [];
JE = [];

for j = 1:num_z
    JT = [JT; jacobian_integrated_intensity_semiinfinite_space_CW(r_1, r_2, r_3, mua, c, mus, g, isRelative)];
    JE = [JE; jacobian_moment_semiinfinite_space_FD(r_1, r_2, r_3, mua, c, mus, g, order, isMellin)];
    r_2 = [r_2(1) 0 j/20];
end


subplot(2,1,1), hold on,
yyaxis left
plot((1:num_z)/20, JT(:,1), 'linewidth', 2)
ylabel('J^{(T)}_\alpha')
yyaxis right
plot((1:num_z)/20, JE(:,1), 'linewidth', 2)
xlabel('z')
ylabel('J^{(\Gamma)}_\alpha')
grid()

subplot(2,1,2), hold on,
mx = max([JT(:,2); JE(:,2)]);
mn = min([JT(:,2); JE(:,2)]);
z_1 = 1/ ( (1 - g) * mus );
yyaxis left
plot((1:num_z)/20, JT(:,2), 'linewidth', 2)
ylabel('J^{(T)}_\nu')
yyaxis right
plot((1:num_z)/20, JE(:,2), 'linewidth', 2)
ylabel('J^{(\Gamma)}_{\nu}')
xlabel('z')
grid()

