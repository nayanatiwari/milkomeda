%% GALAXIES
set (0, 'defaultfigurecolor', [1 1 1])
set(0, 'defaultAxesFontSize', 14)
set(0, 'defaultfigureposition', [0 0 700 350])
clear all; close all;
format compact
%% Galaxy orbits
close all; clear all;
pc = physicsConstants();

Ma = 2e42; % Andromeda Galaxy mass in kg
Mm = 2e42; % will be Milky Way Galaxy mass 1.5e42 kg
e = 25; %eccentricity

kpc = 3.0857e19; %1 kpc in m

xa0 = -30*kpc;
ya0 = -200*kpc;
xm0 = 30*kpc;
ym0 = 200*kpc;

Bx = (Ma*xa0 + Mm*xm0)/(Ma + Mm);
By = (Ma*ya0 + Mm*ym0)/(Ma + Mm);

vxa0 = 0;
vya0 = sqrt(pc.G*Mm*dist(xa0, ya0, Bx, By)*(1 + e)...
        / dist(xa0, ya0, xm0, ym0) ^ 2);
vxm0 = 0;
vym0 = -sqrt(pc.G*Ma*dist(xa0, ya0, Bx, By)*(1 + e)...
        / dist(xa0, ya0, xm0, ym0) ^ 2);

ax1 = @(x1, y1, x2, y2)(-Mm * pc.G * (x1 - x2)) / (dist(x1, y1, x2, y2))^3;
ay1 = @(x1, y1, x2, y2)(-Mm * pc.G * (y1 - y2)) / (dist(x1, y1, x2, y2))^3;
ax2 = @(x1, y1, x2, y2)(-Ma * pc.G * (x2 - x1)) / (dist(x1, y1, x2, y2))^3;
ay2 = @(x1, y1, x2, y2)(-Ma * pc.G * (y2 - y1)) / (dist(x1, y1, x2, y2))^3;
 
ic = [xa0 xm0 ya0 ym0 vxa0 vxm0 vya0 vym0];
tspan = [0 3.15e7*1e10];
derivs = @(t, curr)[curr(5); curr(6); curr(7); curr(8); ...
    ax1(curr(1), curr(3), curr(2), curr(4));...
    ax2(curr(1), curr(3), curr(2), curr(4));...
    ay1(curr(1), curr(3), curr(2), curr(4));...
    ay2(curr(1), curr(3), curr(2), curr(4))];
[t, data] = ode45(derivs, tspan, ic);

xa = data(:,1) ./ kpc;
ya = data(:,3) ./ kpc;
xm = data(:,2) ./ kpc;
ym = data(:,4) ./ kpc;


plot(xa, ya, 'b-');
hold on
plot(xm, ym, 'r-');
hold on
plot(Bx, By, 'k+', 'MarkerSize', 15);

%% Stars
close all; clear all;
M = 2e42;
kpc = 3.0857e19;
e = 0; %explodes for anything other than zero
pc = physicsConstants();
xStars = [];
yStars = [];
vxStars = [];
vyStars = [];
theta = [];
xGal = 150*kpc; %m
yGal = 200*kpc;
vxGal = 100e3; %m/s
vyGal = -200e3;

% x, y initial values
for num=5:5:20
    r = num*kpc;
    theta_spacing = linspace(0, 2*pi, num + 1);
    theta_spacing = theta_spacing(1:(end-1));
    theta = [theta, theta_spacing];
    xStars = [xStars, (cos(theta_spacing) .* r)];
    yStars = [yStars, (sin(theta_spacing) .* r)];
end

xStars = xStars + xGal;
yStars = yStars + yGal;

% vx, vy initial values
vxStars = cos(theta + pi/2) .* (sqrt(pc.G*M*(1+e) ./...
    dist(xStars, yStars, xGal, yGal))) + vxGal;
vyStars = sin(theta + pi/2) .* (sqrt(pc.G*M*(1+e) ./...
    dist(xStars, yStars, xGal, yGal))) + vyGal;

subplot(1, 2, 1);
plot(xStars/kpc, yStars/kpc, 'r*');
hold on
plot(xGal/kpc, yGal/kpc, 'k+');
hold on
quiver(xStars/kpc, yStars/kpc, vxStars/kpc, vyStars/kpc, 'b-');

ax = @(x, y, xGal, yGal) -pc.G*M.*(x - xGal)./...
    dist(x, y, xGal, yGal).^3;
ay = @(x, y, xGal, yGal) -pc.G*M.*(y - yGal)./...
    dist(x, y, xGal, yGal).^3;

tspan = [0 1e9*3.154e7];
ic = [xStars, yStars, vxStars, vyStars, ...
    xGal, yGal, vxGal, vyGal];
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
derivs = @(t,curr) [curr(101:150); curr(151:200); ...
    ax(curr(1:50), curr(51:100), curr(201), curr(202)); ...
    ay(curr(1:50), curr(51:100), curr(201), curr(202)); ...
     curr(203); curr(204); 0; 0];

[t, data] = ode45(derivs, tspan, ic);

xdata = data(:,1:50)/kpc;
ydata = data(:,51:100)/kpc;
vxdata = data(:,101:150)/kpc;
vydata = data(:,151:200)/kpc;

xgdata = data(:,201)/kpc;
ygdata = data(:,202)/kpc;
vxgdata = data(:,203)/kpc;
vygdata = data(:,204)/kpc;

subplot(1, 2, 2);

p1 = plot(xdata(1,:), ydata(1,:), 'r*');
axis([100 500 0 400]);
hold on
p2 = quiver(xdata(1,:), ydata(1,:),...
   vxdata(1,:), vydata(1,:), 'b-');
hold on
p3 = plot(xgdata(1), ygdata(1), 'k+');
hold on
p4 = quiver(xgdata(1), ygdata(1), ...
    vxgdata(1), vygdata(1), 'b-');
hold on

for i=2:length(xdata)
    set(p1, 'XData', xdata(i,:), 'YData', ydata(i,:));
    set(p2, 'XData', xdata(i,:), 'YData', ydata(i,:), ...
        'UData', vxdata(i,:), 'VData', vydata(i,:));
    set(p3, 'XData', xgdata(i), 'YData', ygdata(i));
    set(p4, 'XData', xgdata(i), 'YData', ygdata(i), ...
        'UData', vxgdata(i), 'VData', vygdata(i));
    pause(0.01);
end

%% Stars and Galaxies
close all; clear all;
Msun = 1.989e30; %solar mass in kg
Ma = 1e11 * Msun; % Andromeda Galaxy mass is 2e42 kg
Mm = Ma; % will be Milky Way Galaxy mass 1.5e42 kg
kpc = 3.0857e19;
e = 0; %eccentricity, explodes for anything other than 0
eGal = 30; %eccentricity
pc = physicsConstants();
xStars = [];
yStars = [];
theta = [];
xGal = 20*kpc; %Milky way
yGal = 300*kpc;

xaGal = -20*kpc; %Andromeda
yaGal = -300*kpc;

Bx = (Ma*xaGal + Mm*xGal)/(Ma + Mm);
By = (Ma*yaGal + Mm*yGal)/(Ma + Mm);

vxGal = 0;
vyGal = -sqrt(pc.G*Ma*dist(xaGal, yaGal, Bx, By)*(1 + eGal)...
        / dist(xaGal, yaGal, xGal, yGal) ^ 2);
vxaGal = 0;
vyaGal = sqrt(pc.G*Mm*dist(xaGal, yaGal, Bx, By)*(1 + eGal)...
        / dist(xaGal, yaGal, xGal, yGal) ^ 2);

% x, y initial values
for num=5:5:20
    r = num*kpc;
    theta_spacing = linspace(0, 2*pi, num + 1);
    theta_spacing = theta_spacing(1:(end-1));
    theta = [theta, theta_spacing];
    xStars = [xStars, (cos(theta_spacing) .* r)];
    yStars = [yStars, (sin(theta_spacing) .* r)];
end

%star positions in kpc
xaStars = xStars + xaGal;
yaStars = yStars + yaGal;
xStars = xStars + xGal;
yStars = yStars + yGal;


% vx, vy initial values
vxStars = cos(theta + pi/2) .* (sqrt(pc.G*Mm*(1+e) ./...
    dist(xStars, yStars, xGal, yGal))) + vxGal;
vyStars = sin(theta + pi/2) .* (sqrt(pc.G*Mm*(1+e) ./...
    dist(xStars, yStars, xGal, yGal))) + vyGal;

%reversed orbut for andromeda galaxy - reasoning for behavior 
% in the linked paper in problem statement
vxaStars = cos(theta - pi/2) .* (sqrt(pc.G*Ma*(1+e) ./...
    dist(xaStars, yaStars, xaGal, yaGal))) + vxaGal;
vyaStars = sin(theta - pi/2) .* (sqrt(pc.G*Ma*(1+e) ./...
    dist(xaStars, yaStars, xaGal, yaGal))) + vyaGal;


% accelerations
ax = @(x, y, xg, yg, xa, ya) -pc.G*Mm.*(x - xg)./...
    dist(x, y, xg, yg).^3 + -pc.G*Ma.*(x - xa)./...
    dist(x, y, xa, ya).^3 ;
ay = @(x, y, xg, yg, xa, ya) -pc.G*Mm.*(y - yg)./...
    dist(x, y, xg, yg).^3 + -pc.G*Ma.*(y - ya)./...
    dist(x, y, xa, ya).^3;

axGal = @(x, y, xa, ya)(-Ma * pc.G * (x - xa)) /...
    (dist(x, y, xa, ya))^3;
ayGal = @(x, y, xa, ya)(-Ma * pc.G * (y - ya)) / ...
    (dist(x, y, xa, ya))^3;

axaGal = @(x, y, xa, ya)(-Mm * pc.G * (xa - x)) / ...
    (dist(x, y, xa, ya))^3;
ayaGal = @(x, y, xa, ya)(-Mm * pc.G * (ya - y)) / ...
    (dist(x, y, xa, ya))^3;


tspan = [0 8e9*3.154e7];
ic = [xStars, yStars, xaStars, yaStars, vxStars, vyStars, vxaStars, ...
     vyaStars, xGal, yGal, vxGal, vyGal, xaGal, yaGal, vxaGal, vyaGal];
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
derivs = @(t,curr) [curr(201:250); curr(251:300); ...
    curr(301:350); curr(351:400);...
    ax(curr(1:50), curr(51:100), curr(401), curr(402),...
    curr(405), curr(406)); ...
    ay(curr(1:50), curr(51:100), curr(401), curr(402), ...
    curr(405), curr(406)); ...
    ax(curr(101:150), curr(151:200), curr(401), curr(402),...
    curr(405), curr(406)); ...
    ay(curr(101:150), curr(151:200), curr(401), curr(402), ...
    curr(405), curr(406)); ...    
     curr(403); curr(404); ...
     axGal(curr(401), curr(402), curr(405), curr(406));...
     ayGal(curr(401), curr(402), curr(405), curr(406));...
     curr(407); curr(408);...
     axaGal(curr(401), curr(402), curr(405), curr(406));...
     ayaGal(curr(401), curr(402), curr(405), curr(406))];

[t, data] = ode45(derivs, tspan, ic);

xdata = data(:,1:50)/kpc;
ydata = data(:,51:100)/kpc;
vxdata = data(:,201:250)/kpc;
vydata = data(:,251:300)/kpc;

xadata = data(:,101:150)/kpc;
yadata = data(:,151:200)/kpc;
vxadata = data(:,301:350)/kpc;
vyadata = data(:,351:400)/kpc;

xgdata = data(:,401)/kpc;
ygdata = data(:,402)/kpc;
vxgdata = data(:,403)/kpc;
vygdata = data(:,404)/kpc;

xgadata = data(:,405)/kpc;
ygadata = data(:,406)/kpc;
vxgadata = data(:,407)/kpc;
vygadata = data(:,408)/kpc;

% plot galaxy paths
plot(xgdata, ygdata, 'b-');
hold on
plot(xgadata, ygadata, 'g-'); 
hold on
axis([-300 300 -300 300]);

% animation station
p1 = plot(xdata(1,:), ydata(1,:), 'r*'); %stars
hold on
p2 = quiver(xdata(1,:), ydata(1,:),...
   vxdata(1,:), vydata(1,:), 'b-');
hold on
p3 = plot(xgdata(1), ygdata(1), 'k+'); %milky way
hold on
p4 = quiver(xgdata(1), ygdata(1), ...
    vxgdata(1), vygdata(1), 'c-');
hold on
p5 = plot(xadata(1), yadata(1), 'g*'); %stars
hold on
p6 = quiver(xadata(1), yadata(1), ...
    vxadata(1), vyadata(1), 'm-');
hold on
p7 = plot(xgadata(1), ygadata(1), 'k+'); %andromeda
hold on
p8 = quiver(xgadata(1), ygadata(1), ...
    vxgadata(1), vygadata(1), 'b-');
hold on
plot(Bx, By, 'k+', 'MarkerSize', 15); % barycenter
hold on

for i=2:length(xdata)
    set(p1, 'XData', xdata(i,:), 'YData', ydata(i,:));
    set(p2, 'XData', xdata(i,:), 'YData', ydata(i,:), ...
        'UData', vxdata(i,:), 'VData', vydata(i,:));
    set(p3, 'XData', xgdata(i), 'YData', ygdata(i));
    set(p4, 'XData', xgdata(i), 'YData', ygdata(i), ...
        'UData', vxgdata(i), 'VData', vygdata(i));
    set(p5, 'XData', xadata(i,:), 'YData', yadata(i,:));
    set(p6, 'XData', xadata(i,:), 'YData', yadata(i,:), ...
        'UData', vxadata(i,:), 'VData', vyadata(i,:));
    set(p7, 'XData', xgadata(i), 'YData', ygadata(i));
    set(p8, 'XData', xgadata(i), 'YData', ygadata(i), ...
        'UData', vxgadata(i), 'VData', vygadata(i));
    pause(0.001);
end

%% functions 
function [d] = dist(x1, y1, x2, y2)
d = sqrt((x2 - x1).^2 + (y2 - y1).^2);
end
