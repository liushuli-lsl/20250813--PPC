%% demo_gppf_compare.m
% Original vs Modified Global PPF (G-PPF)
% - Computes rho(t), drho(t) for both designs
% - Plots envelopes +/- rho(t) and shows more pronounced dynamics
% Yi-friendly defaults & clear "knobs" to tune visibility

clear; clc; close all;

%% ===== 0) Global parameters =====
Tp = 2.0;                    % predefined time
a  = 0.03;                   % steady-state accuracy
p  = 1.0;                    % normalized time exponent
t  = linspace(0.01, 1.6*Tp, 1500);    % time axis

% ---- Original design knobs ----
iota   = 0.5;                % gate smoothness (post-Tp)
Sigma0 = 0.15;               % mild widening after Tp
sigma0 = 1.0;                % constant pre-Tp gain

% ---- Modified design knobs (make dynamics more visible) ----
params = struct();
params.alpha   = 1.2;        % kernel scale
params.gamma   = 2.0;        % kernel exponent (steeper early shrink)
params.nu      = 2.0;        % window steepness (larger -> steeper)
params.kappa   = 10.0/(Tp^2);% harder C1 gate after Tp
params.lambda1 = 2*a;        % saturation-triggered widening height

% adaptive sigma (pre-Tp) & Sigma (post-Tp) sensitivities
params.k_u = 1.0; params.k_d = 1.0; params.k_e = 1.0;

% LPF time constants (smaller -> more responsive -> more visible)
params.tau_u = 0.08; params.tau_d = 0.08; params.tau_e = 0.08;

% saturation trigger (example pattern; you can wire in your v(t), d(t), e(t))
params.ud = 1.5;             % saturation bound
params.m1 = 5.0; params.m2 = 0.5;

% bounds for projection
params.sig_min = 0.2; params.sig_max = 6.0; params.Sigma_max = 6.0;

% initialize LPF states (for demo we synthesize signals below)
ru = 0; rd = 0; re = 0; eta = 0;

%% ===== 1) Synthesize demo signals for post-Tp widening =====
% Here we emulate: control output v(t) with a few bursts beyond ud,
% disturbance |d(t)| and |e(t)| that are small but non-zero after Tp.
v   = zeros(size(t));
d   = zeros(size(t));
e   = zeros(size(t));

% make several saturation bursts of v after Tp
bursts = [0.5 0.3; 1.2 0.4; 2.0 0.25]; % [center width] in seconds after Tp
for k = 1:size(bursts,1)
    c = Tp + bursts(k,1);
    w = bursts(k,2);
    v  = v + 3.0*exp( -((t-c).^2) / (2*w^2) );  % exceeds ud for short time
end

% mild disturbances/errors after Tp
d = 0.2 * (t>=Tp) .* (1 + 0.3*sin(1.3*(t-Tp)));
e = 0.15* (t>=Tp) .* (1 + 0.5*cos(0.9*(t-Tp)));

% Delta u = u(v)-g(v) is not explicitly used here (we use |v|-ud proxy)
Delta_u = max(0, abs(v) - params.ud); % proxy for saturation magnitude

%% ===== 2) Compute rho & drho for both designs =====
rhoO  = zeros(size(t));  drhoO = rhoO;
rhoM  = zeros(size(t));  drhoM = rhoM;
upO   = rhoO;            loO   = rhoO;
upM   = rhoM;            loM   = rhoM;

% For modified design, we integrate LPFs for ru, rd, re and eta
dt = t(2)-t(1);

for k = 1:numel(t)
    tk = t(k);

    % ---------- ORIGINAL ----------
    [rhoO(k), drhoO(k)] = gppf_original(tk, a, Tp, p, iota, Sigma0, sigma0);
    upO(k) = +rhoO(k); loO(k) = -rhoO(k);

    % ---------- MODIFIED ----------
    % update LPFs and triggers (discrete forward Euler for demonstration)
    if k==1
        ru_k = ru; rd_k = rd; re_k = re; eta_k = eta;
    else
        ru_k = ru + dt*( -ru/params.tau_u + abs(Delta_u(k)) / params.tau_u );
        rd_k = rd + dt*( -rd/params.tau_d + abs(d(k))       / params.tau_d );
        re_k = re + dt*( -re/params.tau_e + abs(e(k))       / params.tau_e );

        sat_term = (sign(v(k)-params.ud)+1)*(v(k)-params.ud) ...
                 + (sign(v(k)+params.ud)-1)*(v(k)+params.ud);
        eta_k = eta + dt*( -params.m1*eta + params.m2*max(0, sat_term~=0) );
    end
    ru = ru_k; rd = rd_k; re = re_k; eta = eta_k;

    parsNow = params;
    parsNow.r_u = ru; parsNow.r_d = rd; parsNow.r_e = re;
    parsNow.v   = v(k); parsNow.d  = d(k); parsNow.e = e(k);
    parsNow.Delta_u = Delta_u(k);
    parsNow.Tp = Tp; parsNow.a = a; % for defaults()

    [rhoM(k), drhoM(k)] = gppf_modified(tk, a, Tp, p, parsNow);
    upM(k) = +rhoM(k); loM(k) = -rhoM(k);
end

%% ===== 3) PLOTS =====

% --- Figure 1: rho(t) comparison ---
figure('Units','pixels','Position',[50 50 800 360]);
plot(t, rhoO, 'LineWidth',1.6); hold on;
plot(t, rhoM, '--', 'LineWidth',1.6);
xline(Tp,'k--','LineWidth',1);
xlabel('Time t (s)'); ylabel('\rho(t)');
title('Performance Boundary \rho(t): Original vs Modified');
legend('Original','Modified','T_p','Location','northeast');
grid on; box on;

% --- Figure 2: +/- envelopes (bands) ---
figure('Units','pixels','Position',[50 450 800 360]);
plot(t, upO,  'Color',[0.85 0.55 0.1], 'LineWidth',1.3); hold on;
plot(t, loO,  'Color',[0.85 0.55 0.1], 'LineWidth',1.3);
plot(t, upM,  '--', 'Color',[0.1 0.45 0.85], 'LineWidth',1.3);
plot(t, loM,  '--', 'Color',[0.1 0.45 0.85], 'LineWidth',1.3);
xline(Tp,'k--','LineWidth',1);
xlabel('Time t (s)'); ylabel('Envelope \pm\rho(t)');
title('Upper/Lower Envelopes: Original vs Modified');
legend('+Original','-Original','+Modified','-Modified','T_p','Location','northeast');
grid on; box on;

% --- Figure 3: derivative drho(t) (for BLF cancellation insight) ---
figure('Units','pixels','Position',[900 50 800 360]);
plot(t, drhoO, 'LineWidth',1.4); hold on;
plot(t, drhoM, '--', 'LineWidth',1.4);
xline(Tp,'k--','LineWidth',1);
xlabel('Time t (s)'); ylabel('\dot{\rho}(t)');
title('Boundary Variation Rate \dot{\rho}(t): Original vs Modified');
legend('Original','Modified','T_p','Location','northeast');
grid on; box on;

% --- Figure 4: demo signals (v, |d|, |e|) & sigma/Sigma intuition ---
% (Not used directly in original formula; used in modified one)
figure('Units','pixels','Position',[900 450 800 360]);
yyaxis left;
plot(t, v, 'LineWidth',1.1); hold on;
yline(params.ud,'k--','u_d','LineWidth',1);
ylabel('v(t)');
yyaxis right;
plot(t, abs(d), ':', 'LineWidth',1.2);
plot(t, abs(e), '-.', 'LineWidth',1.2);
xlabel('Time t (s)'); grid on; box on;
title('Demo Signals Driving the Modified Envelope');
legend('v(t)','u_d','|d(t)|','|e(t)|','Location','northeast');

disp('Done. Tune params.gamma, params.nu, params.kappa, params.lambda1, tau''s, k''s to adjust visibility.');

%% ===== Local functions =====

function [rho, drho] = gppf_original(t, a, Tp, p, iota, Sigma0, sigma0)
% Original G-PPF:
%   b=(t/Tp)^p, s=1-3b^2+2b^3, phi=-log(b)
%   rho=a + sigma0*phi*s,   t<Tp
%   rho=a*(1+g*Sigma0),     t>=Tp with g=1-exp(-iota*(t-Tp)^2)
if t<=0, t = 1e-6; end
b  = (t/Tp)^p; b = max(min(b,1.0),1e-6);
s  = 1 - 3*b.^2 + 2*b.^3;
sb = -6*b + 6*b.^2;
phi  = -log(b);
phib = -1/b;
dbdt = (p/Tp)*(t/Tp)^(p-1);

if t < Tp
    rho  = a + sigma0 * phi * s;
    drho = sigma0 * (phib*s + phi*sb) * dbdt;
else
    g  = 1 - exp(-iota*(t-Tp)^2);
    gp = 2*iota*(t-Tp)*exp(-iota*(t-Tp)^2);
    rho  = a*(1 + g*Sigma0);
    drho = a*(gp*Sigma0);
end
end

function [rho, drho] = gppf_modified(t, a, Tp, p, params)
% Modified (steeper + saturation-aware) with safe defaults
% Defaults
params = setdef(params,'alpha',1.2);
params = setdef(params,'gamma',2.0);
params = setdef(params,'nu',2.0);
params = setdef(params,'kappa',10/(Tp^2));
params = setdef(params,'lambda1',2*a);
params = setdef(params,'sigma0',1.0);
params = setdef(params,'k_u',1.0); params = setdef(params,'k_d',1.0); params = setdef(params,'k_e',1.0);
params = setdef(params,'tau_u',0.08); params = setdef(params,'tau_d',0.08); params = setdef(params,'tau_e',0.08);
params = setdef(params,'sig_min',0.2); params = setdef(params,'sig_max',6.0); params = setdef(params,'Sigma_max',6.0);
params = setdef(params,'r_u',0); params = setdef(params,'r_d',0); params = setdef(params,'r_e',0);
params = setdef(params,'eta',0); params = setdef(params,'v',0); params = setdef(params,'ud',1.5);
params = setdef(params,'d',0);  params = setdef(params,'e',0);

if t<=0, t = 1e-6; end
b    = (t/Tp)^p; b = max(min(b,1.0),1e-6);
s    = (1-b)^2 * (1 + params.nu*b);
sb   = -2*(1-b)*(1+params.nu*b) + (1-b)^2*params.nu;
phi  = params.alpha * (-log(b))^params.gamma;
phib = params.alpha * params.gamma * (-log(b))^(params.gamma-1) * (-1/b);
dbdt = (p/Tp)*(t/Tp)^(p-1);

% adaptive sigma with projection
sig = params.sigma0 + params.k_u*params.r_u + params.k_d*params.r_d + params.k_e*params.r_e;
sig = min(max(sig, params.sig_min), params.sig_max);

if t < Tp
    rho  = a + sig * phi * s;
    drho = sig * (phib*s + phi*sb) * dbdt;  % (dot{sigma}若需要可从外部传入加上)
else
    gk  = 1 - exp(-params.kappa*(t-Tp)^2);
    gkp = 2*params.kappa*(t-Tp)*exp(-params.kappa*(t-Tp)^2);
    Sigma = params.lambda1 * tanh(params.eta) + params.k_u*params.r_u + params.k_d*params.r_d + params.k_e*params.r_e;
    Sigma = min(max(Sigma, 0), params.Sigma_max);
    rho  = a * (1 + gk*Sigma);
    dSigma = 0;                      % 若外部有 dot{eta}, dot{r_*}，可传入这里叠加
    drho = a * (gkp*Sigma + gk*dSigma);
end
end

function s = setdef(s,field,val)
if ~isfield(s,field) || isempty(s.(field)), s.(field)=val; end
end
