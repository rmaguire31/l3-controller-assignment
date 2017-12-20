function model_system(G)
%MODEL_SYSTEM model system with transfer function
%   Detailed explanation goes here
s = tf('s');

if nargin < 1
    G = (s + 1)/(s*(s^2 + 4*s + 5));
end

spec = struct(...
    'RiseTime',     {1,0.5}, ... seconds
    'SettlingTime', {4,nan}, ... seconds
    'SettlingMin',  {nan,nan}, ...
    'SettlingMax',  {nan,nan}, ...
    'Overshoot',    {5,20}, ... percent 
    'Undershoot',   {nan,nan}, ... percent
    'Peak',         {nan,nan}, ...
    'PeakTime',     {nan,nan}, ... seconds
    'SSError',      {nan,5} ... percent
);

lnM0 = log([spec.Overshoot]./100);
zeta_min = sqrt(lnM0.^2./(pi^2 + lnM0.^2));

zeta = arrayfun(@(min)linspace(min,1,1000)', zeta_min, 'UniformOutput', false);
zeta = [zeta{:}];

omega_Tr = (2.16*zeta + 0.6)./[spec.RiseTime];
omega_Ts = -log(0.02)./(zeta.*[spec.SettlingTime]);
omega = max(omega_Tr, omega_Ts);

x = [-100*zeta(1,:)
     -omega.*zeta
     -flipud(omega.*zeta)
     -100*zeta(1,:)];
y = [100*sqrt(1-zeta(1,:).^2)
     omega.*sqrt(1-zeta.^2)
     -flipud(omega.*sqrt(1-zeta.^2))
     -100*sqrt(1-zeta(1,:).^2)];

k = linspace(0,10,1000);
[r] = rlocus(G,k);
zetar = cos(pi-angle(r(3,:)));
omegar = abs(r(3,:));
k_min = cell(size(spec));
k_max = cell(size(spec));
for i=1:length(spec)
    omegai = interp1(zetar, omegar, zeta(:,i));
    ki = interp1(zetar, k, zeta(:,i));
    k_min{i} = min(ki(omegai>omega(:,i)));
    k_max{i} = max(ki(omegai>omega(:,i)));

    figure();
    rlocus(G);
    hold('on');
    fill(x(:,i),y(:,i), 'c');
    hold('off');
end
disp(k_min);
disp(k_max);

Gc = mean([k_min{:};k_max{:}]');
for i = 1:length(spec)
    sys = feedback(series(Gc(i), G), 1);
    figure();
    step(sys);
    info(i) = stepinfo(sys);
    sys
    disp(info(i))
end
end
