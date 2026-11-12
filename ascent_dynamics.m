function [dx, ceq] = ascent_dynamics(x, u, ~, t, vdat)
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
%  3D Dynamics (rotating sphere) of the ascent phase of a launcher (3S).
%  Optimized version to help convergence:
%   - sin/cos/tan are only computed once as well as the constant terms
%   - T and mdot -> T1 + dT12*s1 + dT23*s2 (less computations)
%   - r >= R_earth is limited by the bounds (xl(1) = R_earth)
%
%  State:   x = [r, theta, phi, v, gamma, chi, m_state]  [km, rad, rad, km/s, rad, rad, t]
%  Control: u = [ux, uy, uz]  (Right-handed UEN unit thrust vector)
%
%  m_state = mass "without jumps": only decreasing for the propellant burnt.
%  Physical mass m = m_state - ms1*s1(t) - ms2*s2(t), with s_i scaling tanh()
% =========================================================================

mu      = vdat.mu;
g0      = vdat.g0;          % [km/s^2]
R_earth = vdat.R_earth;
H       = vdat.H;
omega   = 7.2921e-5;        % [rad/s]

%% --- SCALAR CONSTANTS ---
T1  = vdat.T(1);  T2  = vdat.T(2);  T3  = vdat.T(3);
md1 = T1/(vdat.Isp(1)*g0);
md2 = T2/(vdat.Isp(2)*g0);
md3 = T3/(vdat.Isp(3)*g0);
kD  = 0.5*vdat.rho0*vdat.Cd*vdat.A_ref;      % 0.5*rho0*Cd*A_ref
w2  = omega^2;
tw  = 2*omega;

r     = x(:,1);  phi   = x(:,3);  v   = x(:,4);
gamma = x(:,5);  chi   = x(:,6);  mst = x(:,7);

% Right-handed UEN control components
uU    = u(:,1);  % Up
uE    = u(:,2);  % East
uN    = u(:,3);  % North
t     = t(:) .* ones(size(r));

%% --- STAGING ---
s1 = 0.5*(1 + tanh((t - vdat.t1)/vdat.w));   % 0 -> 1 at separation 1
s2 = 0.5*(1 + tanh((t - vdat.t2)/vdat.w));   % 0 -> 1 at separation 2

T    = T1  + (T2-T1)*s1   + (T3-T2)*s2;
mdot = md1 + (md2-md1)*s1 + (md3-md2)*s2;
m    = mst - vdat.ms(1)*s1 - vdat.ms(2)*s2;  % Physical mass

%% --- TRIGONOMETRIC FUNCTIONS ---
sg = sin(gamma);  cg = cos(gamma);
sc = sin(chi);    cc = cos(chi);
sp = sin(phi);    cp = cos(phi);

% Horizontal thrust projections along velocity (uh) and cross-track (uc, increasing chi)
uh = uN.*cc + uE.*sc;                        % Along horizontal velocity
uc = uE.*cc - uN.*sc;                        % Normal to horizontal velocity (towards +chi)

%% --- ATMOSPHERE AND DRAG ---
v2  = v.^2 + 1e-8;
v_s = sqrt(v2);
D   = (kD .* exp(-(r - R_earth) ./ H)) .* v2;

%% --- THRUST VECTOR PROJECTIONS ---
T_v     = T .* (uU.*sg + uh.*cg);
T_gamma = T .* (uU.*cg - uh.*sg);
T_chi   = T .* uc;

%% --- GRAVITY, CORIOLIS, CENTRIFUGAL ---
g   = mu ./ r.^2;
twv = tw .* v_s;
w2r = (w2 .* r) .* cp;

Cor_gamma  = (twv .* cp) .* sc;
Cor_chi    = twv .* (sp.*cg - cp.*sg.*cc);
Cent_v     = w2r .* (sg.*cp - cg.*sp.*cc);
Cent_gamma = w2r .* (cg.*cp + sg.*sp.*cc);
Cent_chi   = (w2r .* sp) .* sc;

%% --- EQUATIONS OF MOTION ---
vcg    = v_s .* cg;
vr     = v2 ./ r;
dr     = v_s .* sg;
dtheta = (vcg .* sc) ./ (r .* cp);
dphi   = (vcg .* cc) ./ r;
dv     = (T_v - D)./m - g.*sg + Cent_v;
dgamma = (T_gamma./m + (vr - g).*cg + Cor_gamma + Cent_gamma) ./ v_s;
dchi   = (T_chi./m + ((vr.*cg).*cg).*sc.*(sp./cp) + Cor_chi + Cent_chi) ./ vcg;
dm     = -mdot;

dx  = [dr, dtheta, dphi, dv, dgamma, dchi, dm];
ceq = uU.^2 + uE.^2 + uN.^2;
end