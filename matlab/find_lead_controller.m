function [C, K] = find_lead_controller(s1)
%% FIND_LEAD_CONTROLLER Find phase lead controller
%   Detailed explanation goes here

s = tf('s');
G = (s + 1)/(s*(s^2 + 4*s + 5));

z = zero(G);
p = pole(G);

z1 = real(s1);

th = sum(angle(s1-pole(G))) - sum(angle(s1-zero(G))) - angle(s1-z1);
th = mod(th-pi, 2*pi);

p1 = z1 + imag(s1)/tan(th);

C = zpk(z1, p1, 1);
sys = series(C, G);

K = prod(abs(s1-pole(sys)))/prod(abs(s1-zero(sys)));

figure(1);
rlocus(sys);
end