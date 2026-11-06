function [problem, guess] = ascent(m0_tot, stages, Isp, Thrust, g0, lon0, lat0, chi0)
%  ========================================================================
%  AUTHOR  : Gabriele Borghini - 2026
%  ========================================================================
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
%  ICLOCS2 + IPOPT:
%  3 stage ascent.
%   - fixed separation times: t1 = mp1/mdot1, t2 = t1 + mp2/mdot2
%   - mass, thrust and Isp change as functions of t inside the dynamics
%   - tf is free between [t2, t2 + burn_max_stadio3];
%     costo = min tf (with constant thrust it's equivalent to maximize m_f)
%
%  State:   x = [r, theta, phi, v, gamma, chi, m_state]
%  Control: u = [alpha, beta]
%
%  Isp, Thrust = vectors [1x3] ([s], [kN]);  stages(i).mp, stages(i).ms [ton]
% =========================================================================

%% --- COSTANTS ---
mu      = 398600.4418;       % [km^3/s^2]
R_earth = 6378.1;            % [km]
g0_km   = g0/1000;           % [km/s^2]
rho0    = 1.225e9 /1000;     % [kg/km^3]
H_atm   = 8.44;              % [km]
Cd      = 0.5;
A_ref   = pi*(2.5e-3/2)^2;   % [km^2]
inc     = 35.66;             % [deg] target orbit inclination
r_orb   = R_earth + 300;     % [km]

%% --- STAGING ---
mp = [stages(1).mp, stages(2).mp, stages(3).mp];
ms = [stages(1).ms, stages(2).ms, stages(3).ms];
T_kg  = Thrust(:).'/1000;                        % [kg km/s^2]  (1 N = 1e-3 kg km/s^2)
mdot  = Thrust(:).' ./ (Isp(:).'*g0);            % [kg/s]
tb    = mp ./ mdot;                              % burnouut times
t1    = tb(1);
t2    = t1 + tb(2);
tf_max = t2 + tb(3);
w_sw   = 1.5;                                    % [s] step width of tanh

m_end12 = m0_tot - mp(1) - mp(2);                % m_state at end of S2
m_min   = m0_tot - sum(mp) - 1;

%% --- MODEL ---
InternalDynamics = @ascent_dynamics;
SimDynamics      = @ascent_dynamics;

problem.analyticDeriv.gradCost          = [];
problem.analyticDeriv.hessianLagrangian = [];
problem.analyticDeriv.jacConst          = [];
problem.settings = @settings_ascent;

%% --- TIME ---
problem.time.t0_min = 0;      problem.time.t0_max = 0;      guess.t0 = 0;
problem.time.tf_min = t2 + 1; problem.time.tf_max = tf_max;
tf_guess = t2 + 11;           guess.tf = tf_guess;

%% --- PARAMETERS ---
problem.parameters.pl = [];  problem.parameters.pu = [];  guess.parameters = [];

%% --- STATES ---
r0   = R_earth;
v0   = 5e-3;                 % [km/s] non zero to avoid singularity
gam0 = 89.5*pi/180;          % [rad]

problem.states.x0  = [r0, lon0, lat0, v0, gam0, chi0,     m0_tot];
problem.states.x0l = [r0, lon0, lat0, v0, gam0, chi0-0.3, m0_tot];
problem.states.x0u = [r0, lon0, lat0, v0, gam0, chi0+0.3, m0_tot];

%                    r          theta  phi     v    gamma    chi  m
problem.states.xl  = [R_earth,   -pi,  -pi/2,  0,   -pi/2,    0,  m_min];
problem.states.xu  = [r_orb+20,   pi,   pi/2,  10,   pi/2,   pi,  m0_tot];
problem.states.xfl = [r_orb,     -pi,  -pi/2,  0,    0,       0,  m_min];
problem.states.xfu = [r_orb,      pi,   pi/2,  10,   0,      pi,  m_end12];

problem.states.xErrorTol_local    = [1e-4 1e-7 1e-7 1e-6 1e-6 1e-6 1e-3];
problem.states.xErrorTol_integral = [1e-4 1e-7 1e-7 1e-6 1e-6 1e-6 1e-3];
problem.states.xConstraintTol     = [2.0 0.01 0.01 0.05 0.05 0.05 500/1000];

%% --- INPUTS ---
problem.inputs.N   = 0;
problem.inputs.ul  = [-pi/2, -pi];
problem.inputs.uu  = [ pi/2,  pi];
problem.inputs.rl  = [-0.1, -0.1];
problem.inputs.ru  = [ 0.1,  0.1];
problem.inputs.u0l = [pi/2-0.2, chi0-0.3];
problem.inputs.u0u = [pi/2,     chi0+0.3];
problem.inputs.uConstraintTol = [0.1, 0.1];

%% --- SETPOINTS / CONSTRAINTS ---
problem.setpoints.states = [];
problem.setpoints.inputs = [];

problem.constraints.ng_eq    = 0;
problem.constraints.gTol_eq  = [];
problem.constraints.gl       = [];
problem.constraints.gu       = [];
problem.constraints.gTol_neq = [];
% [cos(i); inertial orbital energy = 0 (circular orbit, r and gamma fixed)]
problem.constraints.bl   = [cosd(inc); 0];
problem.constraints.bu   = [cosd(inc); 0];
problem.constraints.bTol = [5e-4; 1e-3];

%% --- DATA ---
problem.data.mu      = mu;
problem.data.g0      = g0_km;
problem.data.R_earth = R_earth;
problem.data.rho0    = rho0;
problem.data.H       = H_atm;
problem.data.Cd      = Cd;
problem.data.A_ref   = A_ref;
problem.data.Isp     = Isp(:).';
problem.data.T       = T_kg;
problem.data.ms      = ms;
problem.data.t1      = t1;
problem.data.t2      = t2;
problem.data.w       = w_sw;

%% --- GUESS (t=0, t1, t2, tf) ---
guess.time = [0; t1; t2; tf_guess];
guess.states(:,1) = [r0;   R_earth+98;  R_earth+297; r_orb];
guess.states(:,2) = [lon0; lon0+0.03;   lon0+0.10;   lon0+0.11];
guess.states(:,3) = [lat0; lat0;        lat0;        lat0];
guess.states(:,4) = [v0;   2.0;         7.0;         7.3];
guess.states(:,5) = [gam0; 0.8;         0.10;        0.0];
guess.states(:,6) = [chi0; chi0;        chi0;        chi0];
guess.states(:,7) = [m0_tot; m0_tot-mp(1); m_end12; m_end12-mdot(3)*11];
guess.inputs(:,1) = [pi/2; 0.8; 0.2; 0.0];
guess.inputs(:,2) = [chi0; chi0; chi0; chi0];

%% --- HANDLES ---
problem.data.InternalDynamics   = InternalDynamics;
problem.data.functionfg         = @fg;
problem.data.plantmodel         = func2str(InternalDynamics);
problem.functions               = {@L, @E, @f, @g, @avrc, @b};
problem.sim.functions           = SimDynamics;
problem.sim.inputX              = [];
problem.sim.inputU              = 1:length(problem.inputs.ul);
problem.functions_unscaled      = {@L_unscaled, @E_unscaled, @f_unscaled, ...
                                   @g_unscaled, @avrc, @b_unscaled};
problem.data.functions_unscaled = problem.functions_unscaled;
problem.data.ng_eq              = problem.constraints.ng_eq;
problem.constraintErrorTol      = [problem.constraints.gTol_eq, problem.constraints.gTol_neq, ...
                                   problem.constraints.gTol_eq, problem.constraints.gTol_neq, ...
                                   problem.states.xConstraintTol, problem.states.xConstraintTol, ...
                                   problem.inputs.uConstraintTol, problem.inputs.uConstraintTol];

%% ================= COSTS / UN-SCALED CONSTRAINTS =================
function stageCost = L_unscaled(~, ~, ~, ~, ~, t, ~)
    stageCost = 0*t;
end

function boundaryCost = E_unscaled(~, ~, ~, ~, ~, ~, tf, ~)
    boundaryCost = tf;     % min tf  (== max final mass)
end

function bc = b_unscaled(~, xf, ~, ~, ~, ~, ~, ~, varargin)
    mu_    = 398600.4418;
    omega  = 7.2921e-5;
    r = xf(1); phi = xf(3); v = xf(4); gam = xf(5); chi = xf(6);
    vU = v*sin(gam);
    vN = v*cos(gam)*cos(chi);
    vE = v*cos(gam)*sin(chi) + omega*r*cos(phi);
    cos_i = cos(phi)*vE / sqrt(vN^2+vE^2);
    e_orb = vU^2 + vN^2 + vE^2 - mu_/r;      % Inertial vs circular velocity
    bc = [cos_i; e_orb];                     % Target orbit

    if ~isempty(varargin) && ~isempty(varargin{1})
        vi = varargin{1};
        if length(vi) == 2
            options = vi{1}; t_segment = vi{2};
            if (strcmp(options.discretization,'hpLGR') || ...
                strcmp(options.discretization,'globalLGR')) && options.adaptseg==1
                if size(t_segment,1) > size(t_segment,2)
                    bc = [bc; diff(t_segment)];
                else
                    bc = [bc, diff(t_segment)'];
                end
            end
        end
    end
end

%% ================= SCALED WRAPPERS (DO NOT MODIFY) =================
function stageCost = L(x, xr, u, ur, p, t, vdat)
    if isfield(vdat,'Xscale')
        x = scale_variables_back(x, vdat.Xscale, vdat.Xshift);
        u = scale_variables_back(u, vdat.Uscale, vdat.Ushift);
        if ~isempty(xr); xr = scale_variables_back(xr, vdat.Xscale, vdat.Xshift); end
        if ~isempty(ur); ur = scale_variables_back(ur, vdat.Uscale, vdat.Ushift); end
    end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        stageCost = 0*t;
    else
        stageCost = L_unscaled(x, xr, u, ur, p, t, vdat);
    end
end

function boundaryCost = E(x0, xf, u0, uf, p, t0, tf, vdat)
    if isfield(vdat,'Xscale')
        x0 = scale_variables_back(x0', vdat.Xscale, vdat.Xshift);
        xf = scale_variables_back(xf', vdat.Xscale, vdat.Xshift);
        u0 = scale_variables_back(u0', vdat.Uscale, vdat.Ushift);
        uf = scale_variables_back(uf', vdat.Uscale, vdat.Ushift);
    end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        boundaryCost = sum(sum(p(:,end-vdat.mode.np*2+1:end)));
    else
        boundaryCost = E_unscaled(x0, xf, u0, uf, p, t0, tf, vdat);
    end
end

function dx = f(x, u, p, t, vdat)
    f_fn = vdat.InternalDynamics;
    if isfield(vdat,'Xscale')
        x = scale_variables_back(x, vdat.Xscale, vdat.Xshift);
        u = scale_variables_back(u, vdat.Uscale, vdat.Ushift);
    end
    dx = f_fn(x, u, p, t, vdat);
    if isfield(vdat,'Xscale')
        dx = scale_variables(dx, vdat.Xscale, 0);
    end
end

function c = g(x, u, p, t, vdat)
    g_fn = vdat.InternalDynamics;
    ng   = nargout(g_fn);
    if isfield(vdat,'Xscale')
        x = scale_variables_back(x, vdat.Xscale, vdat.Xshift);
        u = scale_variables_back(u, vdat.Uscale, vdat.Ushift);
    end
    if ng == 1
        c = [];
    elseif ng == 2
        [~, c] = g_fn(x, u, p, t, vdat);
    else
        [~, ceq, cneq] = g_fn(x, u, p, t, vdat); c = [ceq, cneq];
    end
    if isfield(vdat,'gFilter'); c(:,vdat.gFilter) = []; end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        c = [c-p(:,end-vdat.mode.np*2+1:end-vdat.mode.np), ...
             c+p(:,end-vdat.mode.np+1:end)];
    end
end

function [dx, c] = fg(x, u, p, t, vdat)
    fg_fn = vdat.InternalDynamics;
    ng    = nargout(fg_fn);
    if isfield(vdat,'Xscale')
        x = scale_variables_back(x, vdat.Xscale, vdat.Xshift);
        u = scale_variables_back(u, vdat.Uscale, vdat.Ushift);
    end
    if ng == 1
        c = []; dx = fg_fn(x, u, p, t, vdat);
    elseif ng == 2
        [dx, c] = fg_fn(x, u, p, t, vdat);
    else
        [dx, ceq, cneq] = fg_fn(x, u, p, t, vdat); c = [ceq, cneq];
    end
    if isfield(vdat,'Xscale')
        dx = scale_variables(dx, vdat.Xscale, 0);
    end
    if isfield(vdat,'gFilter'); c(:,vdat.gFilter) = []; end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        c = [c-p(:,end-vdat.mode.np*2+1:end-vdat.mode.np), ...
             c+p(:,end-vdat.mode.np+1:end)];
    end
end

function cr = avrc(x, u, p, t, data)
    [cr] = addRateConstraint(x, u, p, t, data);
end

function bc = b(x0, xf, u0, uf, p, t0, tf, vdat, varargin)
    bc = b_unscaled(x0, xf, u0, uf, p, t0, tf, vdat, varargin);
    if isfield(vdat,'Xscale') && ~isempty(bc)
        x0 = scale_variables_back(x0', vdat.Xscale, vdat.Xshift);
        xf = scale_variables_back(xf', vdat.Xscale, vdat.Xshift);
        u0 = scale_variables_back(u0', vdat.Uscale, vdat.Ushift);
        uf = scale_variables_back(uf', vdat.Uscale, vdat.Ushift);
        bc = b_unscaled(x0, xf, u0, uf, p, t0, tf, vdat, varargin);
    end
end

function dx = f_unscaled(x, u, p, t, vdat)
    Dyn = vdat.InternalDynamics;
    dx  = Dyn(x, u, p, t, vdat);
end

function c = g_unscaled(x, u, p, t, vdat)
    Dyn = vdat.InternalDynamics;
    ng  = nargout(Dyn);
    if ng == 1
        c = [];
    elseif ng == 2
        [~, c] = Dyn(x, u, p, t, vdat);
    else
        [~, ceq, cneq] = Dyn(x, u, p, t, vdat); c = [ceq, cneq];
    end
end

end