clear; clc; close all;

rho = 1.225;          % air density [kg/m^3]
c   = 0.156;          % chord [m]
b   = 1.0;            % span [m]
S   = b * c;          % planform area [m^2]

k_theta = 11.5;       % torsional stiffness [N*m/rad]
k_z     = 300;        % plunge stiffness [N/m]
m       = 2.3;        % mass [kg]
IEA     = 0.0035;     % mass moment of inertia about EA [kg*m^2]

CLalpha = 2*pi;       % lift-curve slope [1/rad]
epsilon = -0.45;      % elastic axis location
gamma_2 = -0.20;      % CG location for Task 2
gamma_3 = 0.0;        % CG location for Task 3

cbar = c/2;           % semi-chord \tilde{c} [m]

%% ============================================================
%  TASK 1: DIVERGENCE SPEED
% ============================================================
% Lecture convention:
% x_AC = cbar/2 (quarter-chord from leading edge)
% x_EA = (1 + epsilon) * cbar

x_AC = cbar/2;
x_EA = (1 + epsilon) * cbar;
e    = x_EA - x_AC;   % positive distance for this geometry

U_div = sqrt((2 * k_theta) / (rho * S * CLalpha * e));

fprintf('Task 1 - Divergence speed: %.2f m/s\n', U_div);

%% ============================================================
%  TASK 2: FLUTTER SPEED AND FREQUENCY (gamma = -0.20)
% ============================================================
x_theta = gamma_2 - epsilon;   % static imbalance parameter

omega_z     = sqrt(k_z / m);
omega_theta = sqrt(k_theta / IEA);

% Find flutter speed from reduced coalescence criterion
flutter_disc_2 = @(U) local_flutter_discriminant( ...
    U, m, IEA, cbar, x_theta, omega_z, omega_theta, rho, epsilon);

U_scan = linspace(0, 80, 4001);
D_scan = arrayfun(flutter_disc_2, U_scan);

idx = find(D_scan(1:end-1) .* D_scan(2:end) <= 0, 1, 'first');
if isempty(idx)
    error('Task 2: No flutter root found in the search interval.');
end

U_f2 = fzero(flutter_disc_2, [U_scan(idx), U_scan(idx+1)]);

% Eigenvalues nu of the first-order state-space system at flutter
A_f2 = local_state_matrix(U_f2, m, IEA, cbar, x_theta, ...
                          omega_z, omega_theta, rho, epsilon);
nu_f2 = eig(A_f2);

% Positive-imaginary eigenvalues correspond to oscillation frequency
nu_pos_2 = nu_f2(imag(nu_f2) > 0);
if isempty(nu_pos_2)
    [~, ii] = max(imag(nu_f2));
    nu_pos_2 = nu_f2(ii);
end

omega_flutter_2 = mean(abs(imag(nu_pos_2)));   % rad/s
f_flutter_2     = omega_flutter_2 / (2*pi);    % Hz

fprintf('Task 2 - Flutter speed: %.2f m/s\n', U_f2);
fprintf('Task 2 - Flutter frequency: %.2f rad/s (%.2f Hz)\n', ...
        omega_flutter_2, f_flutter_2);

% ---------- Task 2 plots: actual eigenvalues nu ----------
U_vec_2 = linspace(0, 1.2 * U_f2, 400);
nu_store_2 = zeros(4, numel(U_vec_2));

for i = 1:numel(U_vec_2)
    A = local_state_matrix(U_vec_2(i), m, IEA, cbar, x_theta, ...
                           omega_z, omega_theta, rho, epsilon);
    nu = eig(A);
    nu_store_2(:, i) = local_sort_eigs(nu);
end

figure;
plot(U_vec_2, real(nu_store_2).', 'LineWidth', 1.2);
xlim([15 30])
xline(U_f2, '--', 'Flutter speed');
xlabel('Airspeed U [m/s]');
ylabel('Real(\nu) [1/s]');
title('Task 2: Real parts of eigenvalues vs Airspeed');
legend('\nu_1','\nu_2','\nu_3','\nu_4','Location','best');
grid on;

figure;
plot(U_vec_2, imag(nu_store_2).', 'LineWidth', 1.2);
xlim([15 30])
xline(U_f2, '--', 'Flutter speed');
xlabel('Airspeed U [m/s]');
ylabel('Imag(\nu) [1/s]');
title('Task 2: Imaginary parts of eigenvalues vs Airspeed');
legend('\nu_1','\nu_2','\nu_3','\nu_4','Location','best');
grid on;

%% ============================================================
%  TASK 3: FLUTTER SPEED AND FREQUENCY (gamma = 0)
% ============================================================
x_theta = gamma_3 - epsilon;   % static imbalance parameter

flutter_disc_3 = @(U) local_flutter_discriminant( ...
    U, m, IEA, cbar, x_theta, omega_z, omega_theta, rho, epsilon);

U_scan = linspace(0, 80, 4001);
D_scan = arrayfun(flutter_disc_3, U_scan);

idx = find(D_scan(1:end-1) .* D_scan(2:end) <= 0, 1, 'first');
if isempty(idx)
    error('Task 3: No flutter root found in the search interval.');
end

U_f3 = fzero(flutter_disc_3, [U_scan(idx), U_scan(idx+1)]);

A_f3 = local_state_matrix(U_f3, m, IEA, cbar, x_theta, ...
                          omega_z, omega_theta, rho, epsilon);
nu_f3 = eig(A_f3);

nu_pos_3 = nu_f3(imag(nu_f3) > 0);
if isempty(nu_pos_3)
    [~, ii] = max(imag(nu_f3));
    nu_pos_3 = nu_f3(ii);
end

omega_flutter_3 = mean(abs(imag(nu_pos_3)));   % rad/s
f_flutter_3     = omega_flutter_3 / (2*pi);    % Hz

fprintf('Task 3 - Flutter speed (gamma = 0): %.2f m/s\n', U_f3);
fprintf('Task 3 - Flutter frequency: %.2f rad/s (%.2f Hz)\n', ...
        omega_flutter_3, f_flutter_3);

% ---------- Task 3 plots: actual eigenvalues nu ----------
U_vec_3 = linspace(0, 1.2 * U_f3, 400);
nu_store_3 = zeros(4, numel(U_vec_3));

for i = 1:numel(U_vec_3)
    A = local_state_matrix(U_vec_3(i), m, IEA, cbar, x_theta, ...
                           omega_z, omega_theta, rho, epsilon);
    nu = eig(A);
    nu_store_3(:, i) = local_sort_eigs(nu);
end

figure;
plot(U_vec_3, real(nu_store_3).', 'LineWidth', 1.2);
xlim([15 30])
xline(U_f3, '--', 'Flutter speed');
xlabel('Airspeed U [m/s]');
ylabel('Real(\nu) [1/s]');
title('Task 3: Real parts of eigenvalues vs Airspeed');
legend('\nu_1','\nu_2','\nu_3','\nu_4','Location','best');
grid on;

figure;
plot(U_vec_3, imag(nu_store_3).', 'LineWidth', 1.2);
xlim([15 30])
xline(U_f3, '--', 'Flutter speed');
xlabel('Airspeed U [m/s]');
ylabel('Imag(\nu) [1/s]');
title('Task 3: Imaginary parts of eigenvalues vs Airspeed');
legend('\nu_1','\nu_2','\nu_3','\nu_4','Location','best');
grid on;

%% ============================================================
%  LOCAL FUNCTIONS
% ============================================================

function d = local_flutter_discriminant(U, m, IEA, cbar, x_theta, ...
                                        omega_z, omega_theta, rho, epsilon)
    % Reduced 2x2 coalescence criterion used only to locate flutter speed
    M = [m*cbar^2,         m*cbar^2*x_theta;
         m*cbar^2*x_theta, IEA];

    K = [m*cbar^2*omega_z^2,                         2*pi*rho*cbar^2*U^2;
         0,     IEA*omega_theta^2 - 2*(0.5+epsilon)*pi*rho*cbar^2*U^2];

    Ared = -M \ K;

    % Coalescence discriminant for the reduced characteristic equation
    d = trace(Ared)^2 - 4*det(Ared);
end

function Astate = local_state_matrix(U, m, IEA, cbar, x_theta, ...
                                      omega_z, omega_theta, rho, epsilon)
    % First-order state-space matrix for plotting actual eigenvalues nu
    M = [m*cbar^2,         m*cbar^2*x_theta;
         m*cbar^2*x_theta, IEA];

    K = [m*cbar^2*omega_z^2,                         2*pi*rho*cbar^2*U^2;
         0,     IEA*omega_theta^2 - 2*(0.5+epsilon)*pi*rho*cbar^2*U^2];

    Astate = [zeros(2), eye(2);
              -M \ K,   zeros(2)];
end

function nu_sorted = local_sort_eigs(nu)
    % Sort by imaginary part first, then by real part, for cleaner plots
    tmp = [imag(nu(:)), real(nu(:))];
    tmp = sortrows(tmp, [-1 2]);
    nu_sorted = complex(tmp(:,2), tmp(:,1));
end
