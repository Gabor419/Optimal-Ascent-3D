%% ========================================================================
%% PROJECT : 3D 3-DOF Optimal Ascent (ICLOCS2 + IPOPT)
%  AUTHOR  : Gabriele Borghini - November 2026
%  Based on ICLOCS2 problem template (Copyright (c) 2018 Imperial College London, MIT License)
% =========================================================================
%  Copyright (c) 2026 Gabriele Borghini. All Rights Reserved.
% 
% COPYRIGHT AND USAGE RIGHTS:
%  The source code, problem formulation, and custom scripts contained herein 
%  are the intellectual property of the author. Permission is hereby granted,
%  free of charge, to any person obtaining a copy of this software and 
%  associated documentation files, to deal in the Software without 
%  restriction, including without limitation the rights to use, copy, 
%  modify, merge, publish, distribute and sublicense copies of the Software,
%  and to permit persons to whom the Software is furnished to do so, subject 
%  to the following conditions:
%  The above copyright notice and this permission notice shall be included 
%  in all copies or substantial portions of the Software.
%
% ACADEMIC CONTEXT & THIRD-PARTY TOOLS:
%  Developed solely for academic research, advanced astrodynamics studies, 
%  and educational purposes at Alma Mater Studiorum - University of Bologna.
%  This implementation interfaces with ICLOCS2 (Imperial College London 
%  Optimal Control Software) and IPOPT, which remain the property of their 
%  respective authors and are subject to their own licensing terms.
%
% DISCLAIMER OF WARRANTY:
%  This software is provided "as is", without warranty of any kind, express
%  or implied, including but not limited to the warranties of accuracy, 
%  numerical stability, or fitness for a particular purpose. In no event 
%  shall the author be liable for any claim or damages arising from the use 
%  of this software.
% -------------------------------------------------------------------------
% DESCRIPTION:
%  This code is an attempt to adapt the problem described in the report
%  for the Launchers and Re-entry course at Alma Mater Studiorum - UniBo,
%  taught by Prof. Marco Sagliano.
% 
%  The original report addressed the problem by solving a TPBVP using primer
%  vector theory (indirect method) via MATLAB's built-in functions.
%  This script instead employs a direct solver (IPOPT via ICLOCS2) to 
%  automatically optimize the trajectory given the constraints and boundary 
%  conditions.
% 
%  Final validation and open-loop verification of the optimal 
%  solution are performed via MATLAB's ode45 integrator with strict tolerances 
%  (RelTol/AbsTol less or equal to 1e-8). 
% 
% NOTE ON GEOGRAPHIC/LAUNCH SITE SELECTION:
%  The choice of the Semnan Space Centre as the launch facility is strictly 
%  driven by physical, geometric, and dynamic problem requirements (such as 
%  latitude, inclination, and trajectory optimization mechanics). This code 
%  is developed solely for academic research and educational purposes, 
%  carrying no political or strategic affiliation.
% =========================================================================

clc; clear; close all; format compact;

if exist('utilities','dir'),    addpath(genpath('utilities'));
elseif exist('../utilities','dir'), addpath(genpath('../utilities')); end
addpath(genpath('D:\OptimalControl\ICLOCS-master\src'))
addpath(genpath('D:\OptimalControl\ICLOCS-master\usr'))
addpath(genpath('D:\OptimalControl\adigator-master'));
delete('Func_*_ascent_dynamics.m');

%% --- VEHICLE DATA ---
m0_tot = 160061/1000;  g0 = 9.80665;
stages(1).mp = 118280/1000; stages(1).ms = 8450.8/1000;   % [t]
stages(2).mp = 25562/1000;  stages(2).ms = 1787.6/1000;
stages(3).mp = 4487.8/1000; stages(3).ms = 337.34/1000;
Isp    = [316, 327, 327];                                 % [s]
Thrust = [1830600/1000, 581040/1000, 613896.29/1000];     % [kN]

% Semnan SC: lat 35.23 N, lon 53.92 E, azimuth 84.08 deg
lat0 = deg2rad(35.23);  lon0 = deg2rad(53.92);  chi0 = deg2rad(84.08);
R_earth = 6378.1;     % [km]

%% --- PROBLEM ---
[problem, guess] = ascent(m0_tot, stages, Isp, Thrust, g0, lon0, lat0, chi0);
options = problem.settings(100);        % 100 points (mesh refinement)
profile on;
[solution, MRHistory] = solveMyProblem(problem, guess, options);
profile off;

%% --- OPEN-LOOP VALIDATION (ode45) ---
vd   = problem.data;
Tn   = solution.T(:);  Un = solution.U;
uRaw = @(t) interp1(Tn, Un, t, 'linear', 'extrap');
uFun = @(t) uRaw(t) ./ sqrt(sum(uRaw(t).^2, 2));
rhs  = @(t,x) ascent_dynamics(x.', uFun(t), [], t, vd).';

odeopt = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,7));
[tv, xv] = ode45(rhs, [Tn(1) Tn(end)], solution.X(1,:).', odeopt);
uv = uFun(tv);

%% --- RESULTS ---
t1 = vd.t1;  t2 = vd.t2;
s1 = @(t) 0.5*(1+tanh((t-t1)/vd.w));
s2 = @(t) 0.5*(1+tanh((t-t2)/vd.w));

% Physical mass subtraction
m_phys_sol = solution.X(:,7) - vd.ms(1)*s1(solution.T) - vd.ms(2)*s2(solution.T);
m_phys_sim = xv(:,7)         - vd.ms(1)*s1(tv)          - vd.ms(2)*s2(tv);

tf = solution.tf;
fprintf('\n=== ASCENT RESULTS ===\n');
fprintf('tf               : %.3f s  (S3 burn: %.3f s)\n', tf, tf - t2);
fprintf('Final altitude   : %.3f km\n', solution.X(end,1) - R_earth);
fprintf('Final latitude   : %.4f deg (Launch: %.4f deg)\n', rad2deg(solution.X(end,3)), rad2deg(lat0));
fprintf('Rel. velocity    : %.4f km/s\n', solution.X(end,4));
fprintf('Final gamma      : %.4f deg\n', rad2deg(solution.X(end,5)));
fprintf('Final azimuth chi: %.4f deg\n', rad2deg(solution.X(end,6)));
fprintf('Final mass (opt) : %.1f kg   [report: 3868.2 kg]\n', m_phys_sol(end)*1000);
fprintf('Final mass (sim) : %.1f kg\n', m_phys_sim(end)*1000);

% Ideal Delta-V 
dV3 = Isp(3)*g0*log((m0_tot - sum([stages(1).mp stages(2).mp]) - sum([stages(1).ms stages(2).ms])) / m_phys_sol(end));
fprintf('dV 3rd stage     : %.3f km/s  [report: 1.253-1.396]\n', dV3/1000);

%% --- DATA REFINEMENT ---
t_opt     = solution.T;
r_opt     = solution.X(:,1);
theta_opt = solution.X(:,2); % Longitude
phi_opt   = solution.X(:,3); % Latitude
v_opt     = solution.X(:,4);
gamma_opt = solution.X(:,5);
chi_opt   = solution.X(:,6);

% Separation events' indexes (computed BEFORE using them for thrust arrays)
[~, idx_sep1] = min(abs(t_opt - t1));
[~, idx_sep2] = min(abs(t_opt - t2));

% Smooth thrust profile (consistent with dynamics)
T_opt = Thrust(1) + (Thrust(2)-Thrust(1))*s1(t_opt) + (Thrust(3)-Thrust(2))*s2(t_opt);

% Velocity components (Local Horizon Frame UEN)
v_up    = v_opt .* sin(gamma_opt);
v_east  = v_opt .* cos(gamma_opt) .* sin(chi_opt);
v_north = v_opt .* cos(gamma_opt) .* cos(chi_opt);

alt_opt = r_opt - R_earth;
alt_sim = xv(:,1) - R_earth;

%% --- KINEMATICS AND MASS PLOTS ---
figure('Name', 'Kinematics overview', 'Color', 'w');

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

% 2. Horizontal Velocities vs Time
subplot(2, 2, 2); hold on;
plot(t_opt, v_east, 'g', 'LineWidth', 2);
plot(t_opt, v_north, 'r', 'LineWidth', 2);
xline(t1, 'k--', 'HandleVisibility', 'off');
xline(t2, 'k--', 'HandleVisibility', 'off');
grid on; xlabel('Time [s]'); ylabel('Velocity [km/s]');
legend('East (v_E)', 'North (v_N)', 'Location', 'best');
title('Horizontal Velocity vs Time'); 
theme 'light';

% 3. Vertical Velocity vs Time
subplot(2, 2, 3); hold on;
plot(t_opt, v_up, 'b', 'LineWidth', 2);
xline(t1, 'k--', 'HandleVisibility', 'off');
xline(t2, 'k--', 'HandleVisibility', 'off');
grid on; xlabel('Time [s]'); ylabel('Velocity [km/s]');
title('Vertical Velocity (v_U) vs Time'); 
theme 'light';

% 4. Pitch Angle vs Time
subplot(2, 2, 4); hold on;
plot(t_opt, rad2deg(asin(max(-1, min(1, solution.U(:,1))))), 'm', 'LineWidth', 2);
xline(t1, 'k--', 'Stage 1 Sep', 'LabelHorizontalAlignment', 'left');
xline(t2, 'k--', 'Stage 2 Sep', 'LabelHorizontalAlignment', 'left');
grid on; xlabel('Time [s]'); ylabel('Pitch \alpha [deg]');
title('Pitch Angle vs Time'); 
theme 'light';

figure('Name', 'Thrust and Mass overview', 'Color', 'w');
subplot(2, 1, 1); hold on;
plot(t_opt, T_opt, 'r', 'LineWidth', 2);
xline(t1, 'k--', 'Stage 1 Sep', 'LabelHorizontalAlignment', 'left');
xline(t2, 'k--', 'Stage 2 Sep', 'LabelHorizontalAlignment', 'left');
grid on; xlabel('Time [s]'); ylabel('Thrust [kN]');
title('Thrust vs Time');
theme 'light';

subplot(2, 1, 2); hold on;
plot(t_opt, m_phys_sol, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
xline(t1, 'k--', 'Stage 1 Sep', 'LabelHorizontalAlignment', 'left');
xline(t2, 'k--', 'Stage 2 Sep', 'LabelHorizontalAlignment', 'left');
grid on; xlabel('Time [s]'); ylabel('Mass [tons]');
title('Vehicle Physical Mass vs Time'); 
theme 'light';

%% --- COSTATES (CARTESIAN), CONTROLS AND HAMILTONIAN ---
if isfield(solution, 'multipliers') && isfield(solution.multipliers, 'lambda_1toN')
    lambda_sph = solution.multipliers.lambda_1toN;
    N_opt = length(t_opt);
    
    lambda_r_cart = zeros(N_opt, 3);
    lambda_v_cart = zeros(N_opt, 3);
    
    % Spherical -> Cartesian (ECEF) via Right-Handed UEN Frame (det = +1)
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
        
        % --- A) Velocity Costates in UEN ---
        cos_gam = cos(gamk);
        if abs(cos_gam) < 1e-6; cos_gam = 1e-6; end
        
        lv_U = lv * sin(gamk) + (lgam / vk) * cos(gamk);
        lv_E = (lv * cos(gamk) - (lgam / vk) * sin(gamk)) * sin(chik) + (lchi / (vk * cos_gam)) * cos(chik);
        lv_N = (lv * cos(gamk) - (lgam / vk) * sin(gamk)) * cos(chik) - (lchi / (vk * cos_gam)) * sin(chik);
        
        % Right-Handed UEN -> ECEF rotation matrix (det(R_uen2ecef) = +1)
        R_uen2ecef = [cos(phik)*cos(thk), -sin(thk), -sin(phik)*cos(thk);
                      cos(phik)*sin(thk),  cos(thk), -sin(phik)*sin(thk);
                      sin(phik),           0,         cos(phik)];
                  
        lambda_v_cart(k, :) = (R_uen2ecef * [lv_U; lv_E; lv_N])';
        
        % --- B) Position Costates in UEN ---
        v_U = vk * sin(gamk);
        v_E = vk * cos(gamk) * sin(chik);
        v_N = vk * cos(gamk) * cos(chik);
        
        cos_phi = cos(phik);
        if abs(cos_phi) < 1e-6; cos_phi = 1e-6; end
        
        lr_U = lr;
        lr_E = (lth / cos_phi - lv_U * v_E + lv_E * v_U + lv_N * v_E * tan(phik) - lv_E * v_N * tan(phik)) / rk;
        lr_N = (lphi - (lv_U * v_N - lv_N * v_U)) / rk;
        
        lambda_r_cart(k, :) = (R_uen2ecef * [lr_U; lr_E; lr_N])';
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
    
    figure('Name', 'Controls vs Time', 'Color', 'w'); hold on;
    u_norm = sqrt(sum(solution.U.^2, 2));
    plot(t_opt, solution.U(:,1), 'b', 'LineWidth', 2);
    plot(t_opt, solution.U(:,2), 'g', 'LineWidth', 2);
    plot(t_opt, solution.U(:,3), 'r', 'LineWidth', 2);
    plot(t_opt, u_norm, 'k:', 'LineWidth', 1.5);
    xline(t1, 'k--', 'Stage 1 Sep', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    xline(t2, 'k--', 'Stage 2 Sep', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    grid on; xlabel('Time [s]'); ylabel('Control components');
    legend('u_U (Up)', 'u_E (East)', 'u_N (North)', '||u||', 'Location', 'best');
    title('Thrust Direction Controls in Right-Handed UEN Frame');
    theme 'light';
    
    % --- C) Hamiltonian ---
    H_tot = zeros(N_opt, 1);
    dyn_func = @ascent_dynamics; 
    n_states = size(lambda_sph, 2); 
    
    for i = 1:N_opt
        dx_i = dyn_func(solution.X(i,:), solution.U(i,:), [], t_opt(i), problem.data);
        dx_i_row = reshape(dx_i(1:n_states), 1, n_states);
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
    fprintf('Multipliers not found in the solution.\n');
end

%% --- 3D TRAJECTORY PLOT (RIGHT-HANDED LOCAL CURVILINEAR ENU) ---
r_sim     = xv(:,1);
theta_sim = xv(:,2); 
phi_sim   = xv(:,3); 

% Curvilinear coordinates wrt the Launch Pad
East  = R_earth .* cos(lat0) .* (theta_opt - lon0);
North = R_earth .* (phi_opt - lat0);
Up    = r_opt - R_earth; 

East_sim  = R_earth .* cos(lat0) .* (theta_sim - lon0);
North_sim = R_earth .* (phi_sim - lat0);
Up_sim    = r_sim - R_earth; 

% Plotted in Right-Handed ENU (X = East, Y = North, Z = Up)
figure('Name', '3D Local Trajectory (East-North-Up)', 'Color', 'w'); hold on;
plot3(East, North, Up, 'r', 'LineWidth', 2.5);
plot3(East_sim, North_sim, Up_sim, '--', 'LineWidth', 2.5, 'Color', [0.5 0.5 0.5]);

plot3(0, 0, 0, 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
plot3(East(idx_sep1), North(idx_sep1), Up(idx_sep1), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
plot3(East(idx_sep2), North(idx_sep2), Up(idx_sep2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
plot3(East(end), North(end), Up(end), 'kp', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

grid on;
xlabel('East [km]'); ylabel('North [km]'); zlabel('Up (Altitude) [km]');
title('3D Local Trajectory (Right-Handed Curvilinear ENU - Semnan Origin)');
legend('Trajectory (IPOPT)', 'Trajectory (ode45)', 'Launch Pad (0,0,0)', 'Stage Separations', '', 'Orbit Insertion', 'Location', 'best');
view(-35, 25);
theme 'light';

%% --- 3D TRAJECTORY PLOT (GLOBAL WITH TARGET ORBIT) ---
% ECEF Cartesian coordinates (IPOPT Optimization)
X = r_opt .* cos(phi_opt) .* cos(theta_opt);
Y = r_opt .* cos(phi_opt) .* sin(theta_opt);
Z = r_opt .* sin(phi_opt);

% ECEF Cartesian coordinates (ode45 Simulation)
X_sim = r_sim .* cos(phi_sim) .* cos(theta_sim);
Y_sim = r_sim .* cos(phi_sim) .* sin(theta_sim);
Z_sim = r_sim .* sin(phi_sim);

% --- 1) Achieved Orbit (from actual final state & INERTIAL velocity) ---
r_f   = r_opt(end);
th_f  = theta_opt(end);
phi_f = phi_opt(end);
v_f   = v_opt(end);
gam_f = gamma_opt(end);
chi_f = chi_opt(end);

omega_E = 7.2921e-5; % [rad/s] Earth rotation rate

% Right-Handed UEN -> ECEF transformation matrix at final point (det = +1)
R_uen2f = [cos(phi_f)*cos(th_f), -sin(th_f), -sin(phi_f)*cos(th_f);
           cos(phi_f)*sin(th_f),  cos(th_f), -sin(phi_f)*sin(th_f);
           sin(phi_f),            0,          cos(phi_f)];

% Inertial velocity in UEN frame (including Earth rotation along East)
v_iner_uen = [v_f * sin(gam_f); ...
              v_f * cos(gam_f) * sin(chi_f) + omega_E * r_f * cos(phi_f); ...
              v_f * cos(gam_f) * cos(chi_f)];

V_ecef_f = R_uen2f * v_iner_uen;
R_ecef_f = [X(end); Y(end); Z(end)];

% Specific angular momentum (h = r x v_inertial) to define the achieved orbital plane
h_vec   = cross(R_ecef_f, V_ecef_f);
n_orbit = h_vec / norm(h_vec);

u_inplane1 = R_ecef_f / norm(R_ecef_f);
u_inplane2 = cross(n_orbit, u_inplane1);

theta_orbit = linspace(0, 2*pi, 300);
Orbit_Achieved = zeros(300, 3);
for k = 1:300
    r_point = r_f * (u_inplane1 * cos(theta_orbit(k)) + u_inplane2 * sin(theta_orbit(k)));
    Orbit_Achieved(k, :) = r_point';
end

% --- 2) Theoretical Target Orbit (i = 35.66 deg, r = R_earth + 300 km) ---
inc_target = deg2rad(35.66);
r_target   = R_earth + 300;

sin_uf = max(-1, min(1, sin(phi_f) / sin(inc_target)));
if cos(chi_f) >= 0
    u_f = asin(sin_uf);          % Ascending branch (v_N >= 0)
else
    u_f = pi - asin(sin_uf);     % Descending branch (v_N < 0)
end

delta_lon    = atan2(sin(u_f) * cos(inc_target), cos(u_f));
Omega_target = th_f - delta_lon;

n_target = [sin(inc_target) * sin(Omega_target); ...
           -sin(inc_target) * cos(Omega_target); ...
            cos(inc_target)];
u_node   = [cos(Omega_target); sin(Omega_target); 0];
u_perp   = cross(n_target, u_node);

Orbit_Target = zeros(300, 3);
for k = 1:300
    r_point = r_target * (u_node * cos(theta_orbit(k)) + u_perp * sin(theta_orbit(k)));
    Orbit_Target(k, :) = r_point';
end

inc_achieved = rad2deg(acos(n_orbit(3)));
fprintf('Target Inclination   : %.4f deg\n', rad2deg(inc_target));
fprintf('Achieved Inclination : %.4f deg (Inertial)\n', inc_achieved);

% --- PLOT ---
figure('Name', '3D Trajectory Global with Target Orbit', 'Color', 'w'); hold on;
[X_earth, Y_earth, Z_earth] = sphere(50);
surf(X_earth * R_earth, Y_earth * R_earth, Z_earth * R_earth, ...
    'EdgeColor', 'none', 'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.3);

plot3(X, Y, Z, 'r', 'LineWidth', 2.5);
plot3(X_sim, Y_sim, Z_sim, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
plot3(Orbit_Achieved(:,1), Orbit_Achieved(:,2), Orbit_Achieved(:,3), 'k-.', 'LineWidth', 1.8);
plot3(Orbit_Target(:,1), Orbit_Target(:,2), Orbit_Target(:,3), 'm:', 'LineWidth', 2.2);

plot3(X(1), Y(1), Z(1), 'k^', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
plot3(X(end), Y(end), Z(end), 'kp', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

grid on; axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('3D Ascent Trajectory & Target Orbit (ECEF)');
legend('Earth', 'Trajectory (IPOPT)', 'Trajectory (ode45)', ...
       'Achieved Orbit (from Insertion)', 'Theoretical Target Orbit (i = 35.66^\circ)', ...
       'Launch Pad', 'Orbit Insertion', 'Location', 'best');
view(60, 20);
theme 'light';