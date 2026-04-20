% ===============================
% TASK 1: Divergence Speed
% ===============================

clear; clc;

% Given parameters
rho = 1.225;          % air density [kg/m^3]
c   = 0.156;          % chord [m]
b   = 1.0;            % span [m]
S   = b * c;          % area [m^2]

k_theta = 11.5;       % torsional stiffness [N*m/rad]
k_z     = 300;        % plunge stiffness [N/m]
m       = 2.3;        % mass [kg]
IEA     = 0.0035;     % mass moment of inertia about EA [kg*m^2]

CLalpha = 2*pi;       % lift slope [1/rad]
epsilon = -0.45;      % elastic axis location
gamma   = -0.20;      % CG location

% Aerodynamic center and elastic axis
x_AC = c/2;
x_EA = (1 + epsilon) * c/2;

% Distance between AC and EA
e = x_AC - x_EA;

% Divergence speed
U_div = sqrt((2 * k_theta) / (rho * S * CLalpha * e));

% Output for Task 1
fprintf('Task 1 - Divergence speed: %.2f m/s\n', U_div);


% ===============================
% TASK 2: Flutter Speed (undamped 2DOF)
% ===============================

cbar = c/2;                 % half-chord
x_theta = gamma - epsilon;  % pitch-plunge coupling parameter

% Uncoupled natural frequencies
omega_z     = sqrt(k_z / m);
omega_theta = sqrt(k_theta / IEA);

% Mass matrix
M = [m*cbar^2,          m*cbar^2*x_theta;
    m*cbar^2*x_theta,  IEA];

% Discriminant of the quadratic eigenvalue problem in lambda = nu^2
% Flutter occurs when the two lambda roots coalesce: discriminant = 0
flutter_disc = @(U) local_flutter_discriminant(U, M, rho, cbar, omega_z, omega_theta, epsilon);

% Find a bracket for the first sign change
U_scan = linspace(0, 80, 4001);
D_scan = arrayfun(flutter_disc, U_scan);
idx = find(D_scan(1:end-1) .* D_scan(2:end) <= 0, 1, 'first');

if isempty(idx)
    error('No flutter root found in the search interval.');
end

% Flutter speed
U_f = fzero(flutter_disc, [U_scan(idx), U_scan(idx+1)]);

% Evaluate eigenvalues at flutter speed
Kf = [m*cbar^2*omega_z^2,                           2*pi*rho*cbar^2*U_f^2;
    0,        IEA*omega_theta^2 - 2*(0.5+epsilon)*pi*rho*cbar^2*U_f^2];

Af = -M \ Kf;
lambda = eig(Af);                       % lambda = nu^2
lambda_flutter = mean(real(lambda));    % repeated root at flutter
omega_flutter = sqrt(-lambda_flutter);  % rad/s
f_flutter = omega_flutter / (2*pi);     % Hz

% Output for Task 2
fprintf('Task 2 - Flutter speed: %.2f m/s\n', U_f);
fprintf('Task 2 - Flutter frequency: %.2f rad/s (%.2f Hz)\n', omega_flutter, f_flutter);


% ===============================
% Local function
% ===============================
function d = local_flutter_discriminant(U, M, rho, cbar, omega_z, omega_theta, epsilon)
K = [M(1,1)*omega_z^2,                         2*pi*rho*cbar^2*U^2;
    0,     M(2,2)*omega_theta^2 - 2*(0.5+epsilon)*pi*rho*cbar^2*U^2];

A = -M \ K;                 % lambda = nu^2
d = trace(A)^2 - 4*det(A);  % discriminant for 2x2 matrix
end

% ===============================
% TASK 2: Plots (Flutter analysis)
% ===============================

U_vec = linspace(0, 1.2*U_f, 400);

lambda1 = zeros(size(U_vec));
lambda2 = zeros(size(U_vec));

nu1 = zeros(size(U_vec));
nu2 = zeros(size(U_vec));

for i = 1:length(U_vec)
    U = U_vec(i);

    % System matrix
    K = [m*cbar^2*omega_z^2,                         2*pi*rho*cbar^2*U^2;
        0,     IEA*omega_theta^2 - 2*(0.5+epsilon)*pi*rho*cbar^2*U^2];

    A = -M \ K;

    lam = eig(A);         % lambda = nu^2

    lambda1(i) = lam(1);
    lambda2(i) = lam(2);

    % convert to nu
    nu = sqrt(lam);

    nu1(i) = nu(1);
    nu2(i) = nu(2);
end

% ===============================
% Plot 1: Eigenvalues lambda
% ===============================
figure;
plot(U_vec, real(lambda1), 'LineWidth', 1.5); hold on;
plot(U_vec, real(lambda2), 'LineWidth', 1.5);
xline(U_f, '--', 'Flutter speed');

xlabel('Airspeed U [m/s]');
ylabel('\lambda = \nu^2');
title('Eigenvalues \lambda vs Airspeed');
legend('\lambda_1','\lambda_2','Location','best');
grid on;


% ===============================
% Plot 2: Frequencies
% ===============================
figure;
plot(U_vec, abs(imag(nu1)), 'LineWidth', 1.5); hold on;
plot(U_vec, abs(imag(nu2)), 'LineWidth', 1.5);
xline(U_f, '--', 'Flutter speed');

xlabel('Airspeed U [m/s]');
ylabel('Frequency [rad/s]');
title('Flutter Frequencies vs Airspeed');
legend('\omega_1','\omega_2','Location','best');
grid on;


% ===============================
% Plot 3: Real part of nu (stability)
% ===============================
figure;
plot(U_vec, real(nu1), 'LineWidth', 1.5); hold on;
plot(U_vec, real(nu2), 'LineWidth', 1.5);
xline(U_f, '--', 'Flutter speed');

xlabel('Airspeed U [m/s]');
ylabel('Real(\nu)');
title('Stability (growth rate) vs Airspeed');
legend('Mode 1','Mode 2','Location','best');
grid on;


% ===============================
% TASK 3: Flutter speed for gamma = 0
% ===============================

gamma = 0;
x_theta = gamma - epsilon;   % dimensionless static imbalance parameter

cbar = c/2;

% Mass matrix
M = [m*cbar^2,        m*cbar^2*x_theta;
     m*cbar^2*x_theta, IEA];

% Uncoupled natural frequencies
omega_z     = sqrt(k_z / m);
omega_theta = sqrt(k_theta / IEA);

% Flutter discriminant
flutter_disc = @(U) local_flutter_discriminant(U, M, rho, cbar, ...
                                               omega_z, omega_theta, epsilon);

% Search for root
U_scan = linspace(0, 80, 4001);
D_scan = arrayfun(flutter_disc, U_scan);
idx = find(D_scan(1:end-1) .* D_scan(2:end) <= 0, 1, 'first');

if isempty(idx)
    error('No flutter root found in the search interval.');
end

U_f3 = fzero(flutter_disc, [U_scan(idx), U_scan(idx+1)]);

% Frequency at flutter
Kf = [m*cbar^2*omega_z^2,                           2*pi*rho*cbar^2*U_f3^2;
      0,        IEA*omega_theta^2 - 2*(0.5+epsilon)*pi*rho*cbar^2*U_f3^2];

Af = -M \ Kf;
lambda = eig(Af);
lambda_flutter = mean(real(lambda));
omega_flutter = sqrt(-lambda_flutter);
f_flutter = omega_flutter / (2*pi);

fprintf('Task 3 - Flutter speed (gamma = 0): %.2f m/s\n', U_f3);
fprintf('Task 3 - Flutter frequency: %.2f rad/s (%.2f Hz)\n', omega_flutter, f_flutter);