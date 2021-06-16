omega = 0; % We want to evaluate moments in the stationary condition
           % i.e. when the full TPSF has been captured
           % i.e. when Lambda(t) when t --> inf
           % which in Fourier space corresponds to Lambda(omega = 0)
           
mua = 0.25;      % [cm^-1]
mus = 20.0;      % [cm^-1]
c   = 0.214e11;  % [cm/s]
g   = 0.72;      % [dimensionless]
isMellin = true; % Only one implemented atm

J = [];

r_1 = [0.0 0.0 0.0]; % [cm] source [x, y, z]
r_2 = [1.5 0.0 0.0]; % [cm] J evaluation [x, y, z]
r_3 = [3.0 0.0 0.0]; % [cm] detector [x, y, z]

for j = 1:100
    J = [J; jacobian_moment_semiinfinite_space_TR(r_1, r_2, r_3, omega, mua, c, mus, g, isMellin)];
    r_2 = [r_2(1) 0 j/20];
end

figure,
subplot(2,1,1), hold on,
mx = max(J(:,1));
mn = min(J(:,1));
z_1 = 1/ ( (1 - g) * mus );
plot((1:100)/20, J(:,1), 'linewidth', 3)
plot([z_1, z_1],[mn, mx], 'linewidth', 2)
legend('J_{\alpha}')
xlabel('z')
ylabel('J(1.5,0,z)')
grid()

subplot(2,1,2), hold on,
mx = max(J(:,2));
mn = min(J(:,2));
z_1 = 1/ ( (1 - g) * mus );
plot((1:100)/20, J(:,2), 'linewidth', 3)
plot([z_1, z_1],[mn, mx], 'linewidth', 2)
legend('J_{\nu}')
xlabel('z')
ylabel('J(1.5,0,z)')
grid()