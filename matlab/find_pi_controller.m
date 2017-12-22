function [Kp,Ki,Kd] = find_pi_controller()
%% FIND_PI_CONTROLLER Find Kp,Ki values of controller
%   Detailed explanation goes here

s = tf('s');
G = (s + 1)/(s*(s^2 + 4*s + 5));

% Perfomance spec
spec = struct(...
    'RiseTime',     {1,0.5}, ... seconds
    'SettlingTime', {4,nan}, ... seconds
    'SettlingMin',  {nan,nan}, ...
    'SettlingMax',  {nan,nan}, ...
    'Overshoot',    {5,20}, ... percent 
    'Undershoot',   {nan,nan}, ... percent
    'Peak',         {nan,nan}, ...
    'PeakTime',     {nan,nan}, ... seconds
    'SSE',          {nan,5} ... percent
);

% Interpret specification
sigmav = -2.5:1e-2:0.5;
omegadv = -10:1e-2:10;
[sigma,omegad] = meshgrid(sigmav,omegadv);
omegan = sqrt(sigma.^2 + omegad.^2);
zeta = abs(sigma)./omegan;
poles = true(length(omegadv),length(sigmav),length(spec));

lnM0 = log([spec.Overshoot]/100);
zeta_M0 = sqrt(lnM0.^2./(lnM0.^2 + pi^2))
poles_M0 = true(size(poles));
for i = 1:length(spec)
    if ~isnan(spec(i).Overshoot)
        poles_M0(:,:,i) = zeta >= zeta_M0(i);
        poles(:,:,i) = poles(:,:,i) & poles_M0(:,:,i);
    end
end

sigma_Tr = log(9)./([spec.RiseTime])
poles_Tr = true(size(poles));
for i = 1:length(spec)
    if ~isnan(spec(i).RiseTime)
        poles_Tr(:,:,i) = sigma <= -sigma_Tr(i);
        poles(:,:,i) = poles(:,:,i) & poles_Tr(:,:,i);
    end
end

sigma_Ts = -log(2/100)./([spec.SettlingTime])
poles_Ts = true(size(poles));
for i = 1:length(spec)
    if ~isnan(spec(i).SettlingTime)
        poles_Ts(:,:,i) = sigma <= -sigma_Ts(i);
        poles(:,:,i) = poles(:,:,i) & poles_Ts(:,:,i);
    end
end

figure(10)
% Root locus
rlocus(G);
hold('on');

% Acceptable region
imagesc(sigmav, omegadv, poles(:,:,1));
colormap([0.9 0.9 0.5;1 1 1]);

% Overshoot condition
xlim = min(sigmav);
ylim = xlim*sqrt(1-zeta_M0(1)^2)/zeta_M0(1);
plot([xlim 0 xlim], [ylim 0 -ylim], '-.k');

% Rise time condition
xlim = [-sigma_Tr(1) -sigma_Tr(1)];
ylim = [min(omegadv) max(omegadv)];
plot(xlim, ylim, '-.k');

% Settling time condition
xlim = [-sigma_Ts(1) -sigma_Ts(1)];
ylim = [min(omegadv) max(omegadv)];
plot(xlim, ylim, '-.k');
hold('off');

% Choose desired zeta and sigma
zeta_crit = zeta_M0;
sigma_crit = max(sigma_Tr, sigma_Ts);

zetap = 1;round(zeta_crit + 0.99*(1 - zeta_crit),2);
omeganp = ceil(sigma_crit./zetap);

% Step response of uncompensated system
[y,t] = step(feedback(G, 1));

% Get time constant and DC gain so we can model it in first order
K = dcgain(feedback(G, 1));
[~,idx] = min(abs(y - K*(1-exp(-1))));
tau = t(idx);

% Plot step response
figure(2);
plot(t,y);
hold('on');
plot([t(1) tau], [y(idx) y(idx)], '-.k');
plot([tau tau], [0 y(idx)], '-.k');
hold('off');

% Detemine controller using zeta/omegan tuning
Kp = (2*zetap.*omeganp*tau-1)/K
Ki = omeganp.^2*tau/K
Kd = zeros(size(spec))
results = cell(size(spec));
targets = cell(size(spec));
for i = 1:length(spec)
    Gc = pid(Kp(i), Kp(i), Kd(i));

    % Analyse step response of compensated system
    sys = feedback(series(Gc, G), 1);
    [y,t] = step(sys);
    figure(2+i);
    plot(t,y);
    info = stepinfo(y, t);
    info.SSE= 100*abs(dcgain(sys) - 1);

    figure(2+length(spec)+i);
    rlocus(series(Gc, G));

    results{i} = struct();
    targets{i} = struct();
    fields = fieldnames(info);
    for j = 1:length(fields)
        if ~isnan(spec(i).(fields{j}))
            results{i}.(fields{j}) = info.(fields{j});
            targets{i}.(fields{j}) = spec(i).(fields{j});
        end
    end
    disp(targets{i});
    disp(results{i});
end
end

