%% ========================================================================
%% PROJECT : 3D 3-DOF Optimal Ascent (ICLOCS2 + IPOPT)
%  AUTHOR  : Gabriele Borghini - 2026
% =========================================================================
%  Copyright (c) 2026 Gabriele Borghini. All Rights Reserved.
%
%  PROPRIETARY AND CONFIDENTIAL:
%  The source code contained herein is the property of Gabriele Borghini.
%  It may not be copied, modified, distributed, or disclosed to any third 
%  party without the express written consent of the author.
%
%  DISCLAIMER:
%  This software is provided "as is", without warranty of any kind, express
%  or implied, including but not limited to the warranties of accuracy or
%  fitness for a particular purpose.
%  ------------------------------------------------------------------------
%  This code is an attempt to adapt the problem described in the report
%  for the Launchers and Re-entry course at Alma Mater Studiorum - UniBo,
%  taught by Prof. Marco Sagliano.
%  The original report addressed the problem by solving a TPBVP using primer
%  vector theory (indirect method) via MATLAB's built-in functions.
%  This script instead employs a direct solver (IPOPT via ICLOCS2) to 
%  automatically optimize the trajectory given the constraints and boundary 
%  conditions.
% 
%  NOTE ON GEOGRAPHIC/LAUNCH SITE SELECTION:
%  The choice of the Semnan Space Centre as the launch facility is strictly 
%  driven by physical, geometric, and dynamic problem requirements (such as 
%  latitude, inclination, and trajectory optimization mechanics). This code 
%  is developed solely for academic research and educational purposes, 
%  carrying no political or strategic affiliation.
%  ========================================================================

clc; clear; close all; format compact;

if exist('utilities','dir'),    addpath(genpath('utilities'));
elseif exist('../utilities','dir'), addpath(genpath('../utilities')); end
addpath(genpath('D:\OptimalControl\ICLOCS-master\src'))
addpath(genpath('D:\OptimalControl\ICLOCS-master\usr'))

%% --- VEHICLE DATA ---
m0_tot = 160061/1000;  g0 = 9.80665;
stages(1).mp = 118280/1000; stages(1).ms = 8450.8/1000;   % [kg/1000]
stages(2).mp = 25562/1000;  stages(2).ms = 1787.6/1000;
stages(3).mp = 4487.8/1000; stages(3).ms = 337.34/1000;
Isp    = [316, 327, 327];                                 % [s]
Thrust = [1830600/1000, 581040/1000, 613896.29/1000];     % [N]

% Semnan SC: lat 35.23 N, lon 53.92 E, azimuth 84.08 deg
lat0 = deg2rad(35.23);  lon0 = deg2rad(53.92);  chi0 = deg2rad(84.08);
R_earth = 6378.1;     % [km]

%% --- PROBLEM ---
[problem, guess] = ascent(m0_tot, stages, Isp, Thrust, g0, lon0, lat0, chi0);
options = problem.settings(600);        % 600 points (mesh refinement)
% options.scaling = 1;                  % Suggested
profile on;
[solution, MRHistory] = solveMyProblem(problem, guess, options);
profile off;
%%
vd   = problem.data;
Tn   = solution.T(:);  Un = solution.U;
uFun = @(t) interp1(Tn, Un, t, 'linear', 'extrap');
rhs  = @(t,x) ascent_dynamics(x.', uFun(t), [], t, vd).';

odeopt = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,7));
[tv, xv] = ode45(rhs, [Tn(1) Tn(end)], solution.X(1,:).', odeopt);
uv = uFun(tv);
%% --- RESULTS ---
vd = problem.data;
t1 = vd.t1;  t2 = vd.t2;
s1 = @(t) 0.5*(1+tanh((t-t1)/vd.w));
s2 = @(t) 0.5*(1+tanh((t-t2)/vd.w));

% phisical mass subtraction
m_phys_sol = solution.X(:,7) - vd.ms(1)*s1(solution.T) - vd.ms(2)*s2(solution.T);
m_phys_sim = xv(:,7)         - vd.ms(1)*s1(tv)          - vd.ms(2)*s2(tv);

tf = solution.tf;
fprintf('\n=== ASCENT RESULTS ===\n');
fprintf('tf               : %.3f s  (burn stadio 3: %.3f s)\n', tf, tf - t2);
fprintf('Final altitude   : %.3f km\n', solution.X(end,1) - R_earth);
fprintf('Rel. velocity    : %.4f km/s\n', solution.X(end,4));
fprintf('Final γ          : %.4f deg\n', rad2deg(solution.X(end,5)));
fprintf('Final mass (opt) : %.1f kg   [report: 3868.2 kg]\n', m_phys_sol(end));
fprintf('Final mass (sim) : %.1f kg\n', m_phys_sim(end));

% Ideal Delta-V 
dV3 = Isp(3)*g0*log((m0_tot - sum([stages(1).mp stages(2).mp]) - sum([stages(1).ms stages(2).ms])) / m_phys_sol(end));
fprintf('ΔV 3rd stage     : %.3f km/s  [report: 1.253-1.396]\n', dV3/1000);

%% --- DATA REFINEMENT ---
% State extraction and derived quantities computation
t_opt     = solution.T;
r_opt     = solution.X(:,1);
theta_opt = solution.X(:,2); % Longitude
phi_opt   = solution.X(:,3); % Latitude
v_opt     = solution.X(:,4);
gamma_opt = solution.X(:,5);
chi_opt   = solution.X(:,6);

% Velocity components (Local Frame)
v_vert  = v_opt .* sin(gamma_opt);
v_north = v_opt .* cos(gamma_opt) .* cos(chi_opt);
v_east  = v_opt .* cos(gamma_opt) .* sin(chi_opt);

alt_opt = r_opt - R_earth;
alt_sim = xv(:,1) - R_earth;

% Separation events' indexes
[~, idx_sep1] = min(abs(t_opt - t1));
[~, idx_sep2] = min(abs(t_opt - t2));

%% --- KINEMATICS AND MASS PLOTS ---
figure('Name', 'Kinematics and Mass overview', 'Color', 'w');

% 1. Altitude vs Time
subplot(2, 2, 1); hold on;
plot(t_opt, alt_opt, 'r', 'LineWidth', 2);
plot(tv, alt_sim, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
xline(t1, 'k--', 'Stage 1 Sep', 'LabelHorizontalAlignment', 'left');
xline(t2, 'k--', 'Stage 2 Sep', 'LabelHorizontalAlignment', 'left');
grid on; xlabel('Time [s]'); ylabel('Altitude [km]');
title('Altitude vs Time'); 
legend('IPOPT', 'ode45', 'Location', 'best');
theme 'light';

% 2. Velocities vs Time
subplot(2, 2, 2); hold on;
plot(t_opt, v_vert, 'b', 'LineWidth', 2); 
plot(t_opt, v_north, 'r', 'LineWidth', 2);
plot(t_opt, v_east, 'g', 'LineWidth', 2);
xline(t1, 'k--', 'HandleVisibility', 'off');
xline(t2, 'k--', 'HandleVisibility', 'off');
grid on; xlabel('Time [s]'); ylabel('Velocity [km/s]');
legend('Vertical', 'North', 'East', 'Location', 'best');
title('Velocity Components vs Time'); 
theme 'light';

% 3. Pitch Angle vs Time
subplot(2, 2, 3); hold on;
plot(t_opt, rad2deg(solution.U(:,1)), 'm', 'LineWidth', 2);
xline(t1, 'k--', 'Stage 1 Sep', 'LabelHorizontalAlignment', 'left');
xline(t2, 'k--', 'Stage 2 Sep', 'LabelHorizontalAlignment', 'left');
grid on; xlabel('Time [s]'); ylabel('Pitch \alpha [deg]');
title('Pitch Angle vs Time'); 
theme 'light';

% 4. Mass vs Time
subplot(2, 2, 4); hold on;
plot(t_opt, m_phys_sol, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
xline(t1, 'k--', 'Stage 1 Sep', 'LabelHorizontalAlignment', 'left');
xline(t2, 'k--', 'Stage 2 Sep', 'LabelHorizontalAlignment', 'left');
grid on; xlabel('Time [s]'); ylabel('Mass [kg]');
title('Vehicle Physical Mass vs Time'); 
theme 'light';

%% --- COSTATES (CARTESIANI) AND HAMILTONIAN ---
if isfield(solution, 'multipliers') && isfield(solution.multipliers, 'lambda_1toN')
    lambda_sph = solution.multipliers.lambda_1toN;
    N_opt = length(t_opt);
    
    lambda_r_cart = zeros(N_opt, 3);
    lambda_v_cart = zeros(N_opt, 3);
    
    % Spherical -> Cartesian (ECEF)
    for k = 1:N_opt
        rk   = r_opt(k);
        thk  = theta_opt(k); 
        phik = phi_opt(k);   
        vk   = max(v_opt(k), 1e-6);
        gamk = gamma_opt(k); 
        chik = chi_opt(k);   
        
        lr   = lambda_sph(k, 1);
        lth  = lambda_sph(k, 2);
        lphi = lambda_sph(k, 3);
        lv   = lambda_sph(k, 4);
        lgam = lambda_sph(k, 5);
        lchi = lambda_sph(k, 6);
        
        % --- A) Velocity Costates ---
        cos_gam = cos(gamk);
        if abs(cos_gam) < 1e-6; cos_gam = 1e-6; end
        
        lv_U = lv * sin(gamk) + (lgam / vk) * cos(gamk);
        lv_N = (lv * cos(gamk) - (lgam / vk) * sin(gamk)) * cos(chik) - (lchi / (vk * cos_gam)) * sin(chik);
        lv_E = (lv * cos(gamk) - (lgam / vk) * sin(gamk)) * sin(chik) + (lchi / (vk * cos_gam)) * cos(chik);
        
        R_loc2ecef = [cos(phik)*cos(thk), -sin(phik)*cos(thk), -sin(thk);
                      cos(phik)*sin(thk), -sin(phik)*sin(thk),  cos(thk);
                      sin(phik),           cos(phik),           0];
                  
        lambda_v_cart(k, :) = (R_loc2ecef * [lv_U; lv_N; lv_E])';
        
        % --- B) Position Costates ---
        v_U = vk * sin(gamk);
        v_N = vk * cos(gamk) * cos(chik);
        v_E = vk * cos(gamk) * sin(chik);
        
        cos_phi = cos(phik);
        if abs(cos_phi) < 1e-6; cos_phi = 1e-6; end
        
        lr_U = lr;
        lr_N = (lphi - (lv_U * v_N - lv_N * v_U)) / rk;
        lr_E = (lth / cos_phi - lv_U * v_E + lv_E * v_U + lv_N * v_E * tan(phik) - lv_E * v_N * tan(phik)) / rk;
        
        lambda_r_cart(k, :) = (R_loc2ecef * [lr_U; lr_N; lr_E])';
    end
    
    figure('Name', 'Cartesian Costates vs Time', 'Color', 'w');
    subplot(2,1,1); hold on;
    plot(t_opt, lambda_r_cart(:,1), 'r', 'LineWidth', 1.5);
    plot(t_opt, lambda_r_cart(:,2), 'g', 'LineWidth', 1.5);
    plot(t_opt, lambda_r_cart(:,3), 'b', 'LineWidth', 1.5);
    xline(t1, 'k--', 'Stage 1', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    xline(t2, 'k--', 'Stage 2', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    grid on; xlabel('Time [s]'); ylabel('\lambda_r [s/km]');
    legend('\lambda_{rx}', '\lambda_{ry}', '\lambda_{rz}', 'Location', 'best');
    title('Position Costates (\lambda_{rx}, \lambda_{ry}, \lambda_{rz})');
    
    subplot(2,1,2); hold on;
    plot(t_opt, lambda_v_cart(:,1), 'r', 'LineWidth', 1.5);
    plot(t_opt, lambda_v_cart(:,2), 'g', 'LineWidth', 1.5);
    plot(t_opt, lambda_v_cart(:,3), 'b', 'LineWidth', 1.5);
    xline(t1, 'k--', 'HandleVisibility', 'off');
    xline(t2, 'k--', 'HandleVisibility', 'off');
    grid on; xlabel('Time [s]'); ylabel('\lambda_v');
    legend('\lambda_{vx}', '\lambda_{vy}', '\lambda_{vz}', 'Location', 'best');
    title('Velocity Costates / Primer Vector (\lambda_{vx}, \lambda_{vy}, \lambda_{vz})');
    theme 'light';
    
    % --- C) Hamiltonian ---
    H_tot = zeros(N_opt, 1);
    
    % Insert here your dynamic function
    dyn_func = @ascent_dynamics; 
    
    % Number of state variables
    n_states = size(lambda_sph, 2); 
    
    for i = 1:N_opt
        % State vector derivative
        dx_i = dyn_func(solution.X(i,:), solution.U(i,:), [], t_opt(i), problem.data);
        
        % Forcing the elements of state vector to be row elements
        dx_i_row = reshape(dx_i(1:n_states), 1, n_states);
        
        % Dot product
        H_tot(i) = dot(lambda_sph(i,:), dx_i_row);
    end
    
    figure('Name', 'Hamiltonian vs Time', 'Color', 'w'); hold on;
    plot(t_opt, H_tot, 'k', 'LineWidth', 2);
    xline(t1, 'k--', 'Stage 1 Sep', 'LabelHorizontalAlignment', 'left');
    xline(t2, 'k--', 'Stage 2 Sep', 'LabelHorizontalAlignment', 'left');
    grid on; xlabel('Time [s]'); ylabel('Hamiltonian H');
    title('Hamiltonian vs Time (Optimization Nodes)'); 
    theme 'light';
    
    h_min = min(H_tot); h_max = max(H_tot);
    if abs(h_max - h_min) < 1e-4
        ylim([h_min - 0.1, h_max + 0.1]);
    else
        ylim([h_min - 0.1*abs(h_min), h_max + 0.1*abs(h_max)]);
    end
else
    fprintf('Moltiplicatori non trovati nella soluzione.\n');
end

%% --- 3D TRAJECTORY PLOT (LOCAL CURVILINEAR NEU) ---
% Spherical coordinate from sim. (ode45)
r_sim     = xv(:,1);
theta_sim = xv(:,2); 
phi_sim   = xv(:,3); 

% Curvilinear coord. wrt the Launch Pad
North = R_earth .* (phi_opt - lat0);
East  = R_earth .* cos(lat0) .* (theta_opt - lon0);
Up    = r_opt - R_earth; 

North_sim = R_earth .* (phi_sim - lat0);
East_sim  = R_earth .* cos(lat0) .* (theta_sim - lon0);
Up_sim    = r_sim - R_earth; 

figure('Name', '3D Local Trajectory (North-East-Up)', 'Color', 'w'); hold on;
plot3(North, East, Up, 'r', 'LineWidth', 2.5);
plot3(North_sim, East_sim, Up_sim, '--', 'LineWidth', 2.5, 'Color', [0.5 0.5 0.5]);

plot3(0, 0, 0, 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
plot3(North(idx_sep1), East(idx_sep1), Up(idx_sep1), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
plot3(North(idx_sep2), East(idx_sep2), Up(idx_sep2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
plot3(North(end), East(end), Up(end), 'kp', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

grid on;
xlabel('North [km]'); ylabel('East [km]'); zlabel('Up (Altitude) [km]');
title('3D Local Trajectory (Curvilinear NEU - Semnan Origin)');
legend('Trajectory (IPOPT)', 'Trajectory (ode45)', 'Launch Pad (0,0,0)', 'Stage Separations', '', 'Orbit Insertion', 'Location', 'best');
view(45, 25);
theme 'light';

%% --- 3D TRAJECTORY PLOT (GLOBALE CON ORBITA TARGET) ---
% Coordinate cartesiane ECEF (Ottimizzazione IPOPT)
X = r_opt .* cos(phi_opt) .* cos(theta_opt);
Y = r_opt .* cos(phi_opt) .* sin(theta_opt);
Z = r_opt .* sin(phi_opt);

% Coordinate cartesiane ECEF (Simulazione ode45)
X_sim = r_sim .* cos(phi_sim) .* cos(theta_sim);
Y_sim = r_sim .* cos(phi_sim) .* sin(theta_sim);
Z_sim = r_sim .* sin(phi_sim);

% --- Generazione Orbita Target passante per il punto finale ---
r_f = r_opt(end);
th_f = theta_opt(end);
phi_f = phi_opt(end);
v_f = v_opt(end);
gam_f = gamma_opt(end);
chi_f = chi_opt(end);

% Matrice di rotazione locale->ECEF al punto finale
R_loc2f = [cos(phi_f)*cos(th_f), -sin(phi_f)*cos(th_f), -sin(th_f);
           cos(phi_f)*sin(th_f), -sin(phi_f)*sin(th_f),  cos(th_f);
           sin(phi_f),           cos(phi_f),           0];
v_local = v_f * [sin(gam_f); cos(gam_f)*cos(chi_f); cos(gam_f)*sin(chi_f)];
V_ecef_f = R_loc2f * v_local;
R_ecef_f = [X(end); Y(end); Z(end)];

% Momento angolare specifico (h = r x v) per definire il piano orbitale
h_vec = cross(R_ecef_f, V_ecef_f);
n_orbit = h_vec / norm(h_vec); % Normale al piano orbitale

% Vettori nel piano orbitale per tracciare l'anello dell'orbita
% Troviamo due vettori ortogonali nel piano dell'orbita
u_inplane1 = R_ecef_f / norm(R_ecef_f);
u_inplane2 = cross(n_orbit, u_inplane1);

% Raggio approssimato dell'orbita target (o semiasse maggiore)
a_target = r_f; 
theta_orbit = linspace(0, 2*pi, 300);
Orbit_ECEF = zeros(300, 3);
for i = 1:300
    % Cerchio/ellisse di raggio r_f nel piano orbitale
    r_point = a_target * (u_inplane1 * cos(theta_orbit(i)) + u_inplane2 * sin(theta_orbit(i)));
    Orbit_ECEF(i, :) = r_point';
end

figure('Name', '3D Trajectory Global with Target Orbit', 'Color', 'w'); hold on;
[X_earth, Y_earth, Z_earth] = sphere(50);
surf(X_earth * R_earth, Y_earth * R_earth, Z_earth * R_earth, ...
    'EdgeColor', 'none', 'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.3);

% Plot delle traiettorie e dell'orbita target
plot3(X, Y, Z, 'r', 'LineWidth', 2.5);
plot3(X_sim, Y_sim, Z_sim, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
plot3(Orbit_ECEF(:,1), Orbit_ECEF(:,2), Orbit_ECEF(:,3), 'k-.', 'LineWidth', 1.5);

% Marker di eventi
plot3(X(1), Y(1), Z(1), 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
% plot3(X(idx_sep1), Y(idx_sep1), Z(idx_sep1), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
% plot3(X(idx_sep2), Y(idx_sep2), Z(idx_sep2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
plot3(X(end), Y(end), Z(end), 'kp', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

grid on; axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Proton-M 3D Ascent Trajectory & Target Orbit (ECEF)');
legend('Earth', 'Trajectory (IPOPT)', 'Trajectory (ode45)', 'Target Orbit', 'Launch Pad', 'Orbit Insertion', 'Location', 'best');
view(60, 20);
theme 'light';