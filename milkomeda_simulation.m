%% Collision of Andromeda and Milky Way Galaxies
% Full simulation in last section
% Nayana Tiwari
set (0, 'defaultfigurecolor', [1 1 1])
set(0, 'defaultAxesFontSize', 14)
set(0, 'defaultfigureposition', [0 0 700 350])
clear all; close all;
format compact
%% Galaxy orbits
close all; clear all;
pc = physicsConstants();
Msun = 1.989e30; %solar mass in kg
Ma = 68e10 * Msun; % Andromeda Galaxy mass is 2e42 kg
Mm = 78e10 * Msun; % will be Milky Way Galaxy mass 1.5e42 kg

kpc = 3.0857e19; %1 kpc in m

xGal = 0; %Milky way
yGal = 0;
zGal = 0;
vxGal = 0;
vyGal = 0;
vzGal = 0;

%andromeda galaxy position
xaGal = -379.2*kpc;
yaGal = 612.7*kpc;
zaGal = -283.1*kpc;

%barycenter
Bx = (Ma*xaGal + Mm*xGal)/(Ma + Mm);
By = (Ma*yaGal + Mm*yGal)/(Ma + Mm);
Bz = (Ma*zaGal + Mm*zGal)/(Ma + Mm);

%center origin on Barycenter
xGal = xGal - Bx;
yGal = yGal - By;
zGal = zGal - Bz;

xaGal = xaGal - Bx;
yaGal = yaGal - By;
zaGal = zaGal - Bz;

%Barycenter on origin
Bx = 0;
By = 0;
Bz = 0;

% from paper
vr = 117e3; %m/s
vt = 42e3; %m/s
va = sqrt(vr^2 + vt^2); %calcualtion check

%angles
phi = atan2(yaGal, xaGal);
ra = sqrt(xaGal^2 + yaGal^2);
theta = asin(zaGal / ra);

%velocity components for andromeda galaxy
vxaGal = vr * sin(theta) * cos(phi) - vt * sin(theta + pi/2) * cos(phi)
vyaGal = vr * sin(theta) * sin(phi) - vt * sin(theta + pi/2) * sin(phi);
vzaGal = vr * cos(theta) - vt * cos(theta + pi/2);
va = sqrt(vxaGal^2 + vyaGal^2 + vzaGal^2);

% transform velocities to immobilize the barycenter
% momentum should be 0
vtot_x = (Ma/(Ma + Mm)) * vxaGal
vtot_y = (Ma/(Ma + Mm)) * vyaGal
vtot_z = (Ma/(Ma + Mm)) * vzaGal

vxaGal = vxaGal - vtot_x
vyaGal = vyaGal - vtot_y;
vzaGal = vzaGal - vtot_z;

vxGal = -1 * vtot_x;
vyGal = -1 * vtot_y;
vzGal = -1 * vtot_z;

% ptot is relatively 0 - yay!!!
ptot_x = Ma*vxaGal + Mm*vxGal
ptot_y = Ma*vyaGal + Mm*vyGal
ptot_z = Ma*vzaGal + Mm*vzGal

ptot = Ma * sqrt(vxaGal^2 + vyaGal^2 + vzaGal^2) - ...
   Mm * sqrt(vxGal^2 + vyGal^2 + vzGal^2)

% ode set up
axGal = @(x1, y1, z1, x2, y2, z2, M)(-M * pc.G * (x1 - x2)) / ...
    (dist(x1, y1, z1, x2, y2, z2))^3;
ayGal = @(x1, y1, z1, x2, y2, z2, M)(-M * pc.G * (y1 - y2)) / ...
    (dist(x1, y1, z1, x2, y2, z2))^3;
azGal = @(x1, y1, z1, x2, y2, z2, M)(-M * pc.G * (z1 - z2)) / ...
    (dist(x1, y1, z1, x2, y2, z2))^3;

ic = [xGal yGal zGal xaGal yaGal zaGal ...
    vxGal vyGal vzGal vxaGal vyaGal vzaGal];
tspan = [0 3.15e7*1e10];
derivs = @(t, curr)[curr(7); curr(8); curr(9); curr(10); ...
    curr(11); curr(12); ...
    % milky way
    axGal(curr(1), curr(2), curr(3), curr(4), curr(5), curr(6), Ma);...
    ayGal(curr(1), curr(2), curr(3), curr(4), curr(5), curr(6), Ma);...
    azGal(curr(1), curr(2), curr(3), curr(4), curr(5), curr(6), Ma);...
    % andromeda galaxy
    axGal(curr(4), curr(5), curr(6), curr(1), curr(2), curr(3), Mm);...
    ayGal(curr(4), curr(5), curr(6), curr(1), curr(2), curr(3), Mm);...
    azGal(curr(4), curr(5), curr(6), curr(1), curr(2), curr(3), Mm)];
[t, data] = ode45(derivs, tspan, ic);

xmw = data(:,1) ./ kpc;
ymw = data(:,2) ./ kpc;
zmw = data(:,3) ./ kpc;
xa = data(:,4) ./ kpc;
ya = data(:,5) ./ kpc;
za = data(:,6) ./ kpc;
vxmw = data(:,7) ./ kpc;
vymw = data(:,8) ./ kpc;
vzmw = data(:,9) ./ kpc;
vxa = data(:,10) ./ kpc;
vya = data(:,11) ./ kpc;
vza = data(:,12) ./ kpc;

figure(100);
plot3(xa(1), ya(1), za(1), 'b+', 'MarkerSize', 15);
hold on
% quiver(xa(1), ya(1), vxa(1) * kpc, sqrt(vya(1)^2 + vza(1)^2)*kpc , 'b-');
plot3(xmw(1), ymw(1),  zmw(1), 'r+', 'MarkerSize', 15);
% quiver(xmw(1), ymw(1), vxmw(1) * kpc, ...
%     sqrt(vymw(1)^2 + vzmw(1)^2)*kpc , 'r-');
hold on
plot3(xa, ya, za, 'b-');
hold on
plot3(xmw, ymw, zmw, 'r-');
hold on
plot3(Bx/kpc, By/kpc, Bz/kpc, 'k+', 'MarkerSize', 15);

%% Stars
close all; clear all;
M = 2e42;
kpc = 3.0857e19;
tilt = pi/3;
e = 0; %explodes for anything other than zero
pc = physicsConstants();

xStars = [];
yStars = [];
zStars = [];
vxStars = [];
vyStars = [];
theta = [];
xGal = 150*kpc; %m
yGal = 200*kpc;
zGal = 0;
vxGal = 100e3; %m/s
vyGal = -200e3;
vzGal = 0;

% x, y initial values
for num=5:5:20
    r = num*kpc;
    theta_spacing = linspace(0, 2*pi, num + 1);
    theta_spacing = theta_spacing(1:(end-1));
    theta = [theta, theta_spacing];

    xStars = [xStars, (cos(theta_spacing).* r)];
    yStars = [yStars, (sin(theta_spacing) .* r)];
end

zStars = zeros(size(xStars));

[xStars, yStars, zStars] = rotationEuler(xStars, yStars, zStars, tilt, 0);

xStars = xStars + xGal;
yStars = yStars + yGal;
zStars = zStars + zGal;

% vx, vy initial values
vxStars = cos(theta + pi/2) .* (sqrt(pc.G*M*(1+e) ./...
    dist(xStars, yStars, zStars, xGal, yGal, zGal))) + vxGal;
vyStars = sin(theta + pi/2) .* (sqrt(pc.G*M*(1+e) ./...
    dist(xStars, yStars, zStars, xGal, yGal, zGal))) + vyGal;
vzStars = zeros(size(vxStars));
[vxStars, vyStars, vzStars] = rotationEuler(vxStars, vyStars,...
   vzStars, tilt, 0);

subplot(1, 2, 1);
plot3(xStars/kpc, yStars/kpc, zStars/kpc, 'r*');
hold on
plot3(xGal/kpc, yGal/kpc, zGal/kpc, 'k+');
hold on
quiver3(xStars/kpc, yStars/kpc, zStars/kpc, vxStars/kpc, vyStars/kpc, ...
    vzStars/kpc, 'b-');
view([45 45]);
 
% phi = atan2(yA, xA);
% rA = sqrt(xA^2 + yA^2);
% theta = asin(zA / rA);

ax = @(x, y, z, xGal, yGal, zGal) -pc.G * M .*(x - xGal) ./...
    dist(x, y, z, xGal, yGal, zGal).^3;
ay = @(x, y, z, xGal, yGal, zGal) -pc.G * M .*(y - yGal) ./...
    dist(x, y, z, xGal, yGal, zGal).^3;
az = @(x, y, z, xGal, yGal, zGal) -pc.G * M .*(z - zGal) ./...
    dist(x, y, z, xGal, yGal, zGal).^3;

tspan = [0 1e9*3.154e7];
ic = [xStars, yStars, zStars, vxStars, vyStars, vzStars, ...
    xGal, yGal, zGal, vxGal, vyGal, vzGal];
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
% x 1:50, y 51:100, z 101:150, vx 151:200, vy 201:250, vz 251:300
% Galaxy: x 301 y 302 z 303 vx 304 vy 305 vz 306
derivs = @(t,curr) [curr(151:200); curr(201:250); curr(251:300);...
    ax(curr(1:50), curr(51:100), curr(101:150), ...
        curr(301), curr(302), curr(303)); ...
    ay(curr(1:50), curr(51:100), curr(101:150),...
        curr(301), curr(302), curr(303)); ...
    az(curr(1:50), curr(51:100), curr(101:150),...
        curr(301), curr(302), curr(303));
     curr(304); curr(305); curr(306); 0; 0; 0];

[t, data] = ode45(derivs, tspan, ic);

xdata = data(:,1:50)/kpc;
ydata = data(:,51:100)/kpc;
zdata = data(:,101:150)/kpc;

vxdata = data(:,151:200)/kpc;
vydata = data(:,201:250)/kpc;
vzdata = data(:,251:300)/kpc;

xgdata = data(:,301)/kpc;
ygdata = data(:,302)/kpc;
zgdata = data(:,303)/kpc;
vxgdata = data(:,304)/kpc;
vygdata = data(:,305)/kpc;
vzgdata = data(:,306)/kpc;

subplot(1, 2, 2);

p1 = plot3(xdata(1,:), ydata(1,:), zdata(1,:),'r*');
hold on
p2 = quiver3(xdata(1,:), ydata(1,:), zdata(1,:), ...
   vxdata(1,:), vydata(1,:), vzdata(1,:), 'b-');
hold on
p3 = plot3(xgdata(1), ygdata(1), zgdata(1), 'k+');
hold on
p4 = quiver3(xgdata(1), ygdata(1), zgdata(1), ...
    vxgdata(1), vygdata(1), vzgdata(1), 'b-');
view([45 45]);
hold on


for i=2:length(xdata)
    % stars
    set(p1, 'XData', xdata(i,:), 'YData', ydata(i,:), 'ZData', zdata(i,:));
    set(p2, 'XData', xdata(i,:), 'YData', ydata(i,:), ... 
        'ZData', zdata(i,:), 'UData', vxdata(i,:), ....
        'VData', vydata(i,:), 'WData', vzdata(i,:));
    % galaxy
    set(p3, 'XData', xgdata(i), 'YData', ygdata(i), 'ZData', zgdata(i));
    set(p4, 'XData', xgdata(i,:), 'YData', ygdata(i,:), ... 
        'ZData', zgdata(i,:), 'UData', vxgdata(i,:), ....
        'VData', vygdata(i,:), 'WData', vzgdata(i,:));
    view([45 45]);
    pause(0.01);
end

%% Full Collision - Stars and Galaxies
close all; clear all;
Msun = 1.989e30; %solar mass in kg
Ma = 68e10 * Msun; % Andromeda Galaxy mass is 2e42 kg
Mm = 78e10 * Msun; % will be Milky Way Galaxy mass 1.5e42 kg
kpc = 3.0857e19;
e = 0; %eccentricity for stars, explodes for anything other than 0
tilt = deg2rad(39.8); %andromeda's tilt from the paper
pc = physicsConstants();

xStars = [];
yStars = [];
zStars = [];
thetaStars = [];

xGal = 0; %Milky way
yGal = 0;
zGal = 0;

%andromeda galaxy position
xaGal = -379.2*kpc;
yaGal = 612.7*kpc;
zaGal = -283.1*kpc;

%barycenter
Bx = (Ma*xaGal + Mm*xGal)/(Ma + Mm);
By = (Ma*yaGal + Mm*yGal)/(Ma + Mm);
Bz = (Ma*zaGal + Mm*zGal)/(Ma + Mm);

%center origin on Barycenter
xGal = xGal - Bx;
yGal = yGal - By;
zGal = zGal - Bz;

xaGal = xaGal - Bx;
yaGal = yaGal - By;
zaGal = zaGal - Bz;

%Barycenter on origin
Bx = 0;
By = 0;
Bz = 0;

% from paper
vr = 117e3; %m/s
vt = 42e3; %m/s
% va = sqrt(vr^2 + vt^2); %calcualtion check

%angles
phi = atan2(yaGal, xaGal);
ra = sqrt(xaGal^2 + yaGal^2);
theta = asin(zaGal / ra);

%velocity components for andromeda galaxy
vxaGal = vr * sin(theta) * cos(phi) - vt * sin(theta + pi/2) * cos(phi)
vyaGal = vr * sin(theta) * sin(phi) - vt * sin(theta + pi/2) * sin(phi);
vzaGal = vr * cos(theta) - vt * cos(theta + pi/2);
va = sqrt(vxaGal^2 + vyaGal^2 + vzaGal^2);

% transform velocities to immobilize the barycenter
% momentum should be 0
vtot_x = (Ma/(Ma + Mm)) * vxaGal;
vtot_y = (Ma/(Ma + Mm)) * vyaGal;
vtot_z = (Ma/(Ma + Mm)) * vzaGal;

vxaGal = vxaGal - vtot_x
vyaGal = vyaGal - vtot_y;
vzaGal = vzaGal - vtot_z;

vxGal = -1 * vtot_x;
vyGal = -1 * vtot_y;
vzGal = -1 * vtot_z;

% stars
for num=5:5:20
    r = 5*num*kpc;
    theta_spacing = linspace(0, 2*pi, num + 1);
    theta_spacing = theta_spacing(1:(end-1));
    thetaStars = [thetaStars, theta_spacing];

    xStars = [xStars, (cos(theta_spacing).* r)];
    yStars = [yStars, (sin(theta_spacing) .* r)];
end

zStars = zeros(size(xStars));

% tilt and set andromeda stars x, y, z
[xaStars, yaStars, zaStars] = rotationEuler(xStars, yStars, zStars, tilt, 0);

%adjust stars to galaxy initial positions
xaStars = xaStars + xaGal;
yaStars = yaStars + yaGal;
zaStars = zaStars + zaGal;

xStars = xStars + xGal;
yStars = yStars + yGal;
zStars = zStars + zGal;

% milky way stars vx, vy initial values
vxStars = cos(thetaStars + pi/2) .* (sqrt(pc.G*Mm*(1+e) ./...
    dist(xStars, yStars, zStars, xGal, yGal, zGal))) + vxGal;
vyStars = sin(thetaStars + pi/2) .* (sqrt(pc.G*Mm*(1+e) ./...
    dist(xStars, yStars, zStars, xGal, yGal, zGal))) + vyGal;
vzStars = zeros(size(vxStars)) + vzGal;


%reversed orbit for andromeda galaxy - reasoning for behavior 
% in the linked paper in problem statement
vxaStars = cos(thetaStars - pi/2) .* (sqrt(pc.G*Ma*(1+e) ./...
    dist(xaStars, yaStars, zaStars, xaGal, yaGal, zaGal)));
vyaStars = sin(thetaStars - pi/2) .* (sqrt(pc.G*Ma*(1+e) ./...
    dist(xaStars, yaStars, zaStars, xaGal, yaGal, zaGal)));
vzaStars = zeros(size(vxaStars));
[vxaStars, vyaStars, vzaStars] = rotationEuler(vxaStars, vyaStars,...
   vzaStars, tilt, 0);

vxaStars = vxaStars + vxaGal;
vyaStars = vyaStars + vyaGal;
vzaStars = vzaStars + vzaGal;


% star accelerations
ax = @(x, y, z, xGal, yGal, zGal, xaGal, yaGal, zaGal)...
    -pc.G * Mm .*(x - xGal) ./ dist(x, y, z, xGal, yGal, zGal).^3 + ...
    -pc.G * Ma .*(x - xaGal) ./ dist(x, y, z, xaGal, yaGal, zaGal).^3;
ay = @(x, y, z, xGal, yGal, zGal, xaGal, yaGal, zaGal)...
    -pc.G * Mm .*(y - yGal) ./ dist(x, y, z, xGal, yGal, zGal).^3 + ...
    -pc.G * Ma .*(y - yaGal) ./ dist(x, y, z, xaGal, yaGal, zaGal).^3;
az = @(x, y, z, xGal, yGal, zGal, xaGal, yaGal, zaGal)...
    -pc.G * Mm .*(z - zGal) ./ dist(x, y, z, xGal, yGal, zGal).^3 + ...
    -pc.G * Ma .*(z - zaGal) ./ dist(x, y, z, xaGal, yaGal, zaGal).^3;

% galaxy accelerations
axGal = @(x1, y1, z1, x2, y2, z2, M)(-M * pc.G * (x1 - x2)) / ...
    (dist(x1, y1, z1, x2, y2, z2))^3;
ayGal = @(x1, y1, z1, x2, y2, z2, M)(-M * pc.G * (y1 - y2)) / ...
    (dist(x1, y1, z1, x2, y2, z2))^3;
azGal = @(x1, y1, z1, x2, y2, z2, M)(-M * pc.G * (z1 - z2)) / ...
    (dist(x1, y1, z1, x2, y2, z2))^3;

tspan = [0 3.15e7*10e11];
% milky way
%     stars: x 1:50,     y 51:100,   z 101:150
%           vx 301:350, vy 351:400, vz 401:450 
%    galaxy: x 601       y 602       z 603
%           vx 607      vy 608      vy 609

% andromeda
%     stars: x 151:200,  y 201:250,  z 251:300
%           vx 451:500, vy 501:550, vz 551:600
%    galaxy: x 604       y 605       z 606
%           vx 610      vy 611      vz 612


ic = [xStars, yStars, zStars, xaStars, yaStars, zaStars, ...
    vxStars, vyStars, vzStars, vxaStars, vyaStars, vzaStars, ...
    xGal, yGal, zGal, xaGal, yaGal, zaGal, ...
    vxGal, vyGal, vzGal, vxaGal, vyaGal, vzaGal];
opts = odeset('RelTol',1e-2,'AbsTol',1e-2);
derivs = @(t,curr) [curr(301:350); curr(351:400); curr(401:450);...
    curr(451:500); curr(501:550); curr(551:600);...
    %milky way stars velocity
    ax(curr(1:50), curr(51:100), curr(101:150), ...
      curr(601), curr(602), curr(603), curr(604), curr(605), curr(606));...
    ay(curr(1:50), curr(51:100), curr(101:150), ...
      curr(601), curr(602), curr(603), curr(604), curr(605), curr(606));...
    az(curr(1:50), curr(51:100), curr(101:150), ...
      curr(601), curr(602), curr(603), curr(604), curr(605), curr(606));...
    %andromeda stars velocity
    ax(curr(151:200), curr(201:250), curr(251:300), ...
      curr(601), curr(602), curr(603), curr(604), curr(605), curr(606));...
    ay(curr(151:200), curr(201:250), curr(251:300), ...
      curr(601), curr(602), curr(603), curr(604), curr(605), curr(606));...
    az(curr(151:200), curr(201:250), curr(251:300), ...
      curr(601), curr(602), curr(603), curr(604), curr(605), curr(606));...
    %galaxy position derivatives
    curr(607); curr(608); curr(609); curr(610); curr(611); curr(612); ...
    % milky way
    axGal(curr(601), curr(602), curr(603), ...
        curr(604), curr(605), curr(606), Ma);...
    ayGal(curr(601), curr(602), curr(603), ...
        curr(604), curr(605), curr(606), Ma);...
    azGal(curr(601), curr(602), curr(603), ...
        curr(604), curr(605), curr(606), Ma);...
    % andromeda galaxy
    axGal(curr(604), curr(605), curr(606), ...
        curr(601), curr(602), curr(603), Mm);...
    ayGal(curr(604), curr(605), curr(606), ...
        curr(601), curr(602), curr(603), Mm);...
    azGal(curr(604), curr(605), curr(606), ...
        curr(601), curr(602), curr(603), Mm)];

[t, data] = ode23(derivs, tspan, ic, opts);

%milky way stars
xdata = data(:,1:50)/kpc;
ydata = data(:,51:100)/kpc;
zdata = data(:,101:150)/kpc;
vxdata = data(:,301:350)/kpc;
vydata = data(:,351:400)/kpc;
vzdata = data(:,401:450)/kpc;

%andromeda stars
xadata = data(:,151:200)/kpc;
yadata = data(:,201:250)/kpc;
zadata = data(:,251:300)/kpc;
vxadata = data(:,451:500)/kpc;
vyadata = data(:,501:550)/kpc;
vzadata = data(:,551:600)/kpc;

%milky way galaxy
xgdata = data(:,601)/kpc;
ygdata = data(:,602)/kpc;
zgdata = data(:,603)/kpc;
vxgdata = data(:,607)/kpc;
vygdata = data(:,608)/kpc;
vzgdata = data(:,609)/kpc;

%andromeda galaxy
xgadata = data(:,604)/kpc;
ygadata = data(:,605)/kpc;
zgadata = data(:,606)/kpc;
vxgadata = data(:,610)/kpc;
vygadata = data(:,611)/kpc;
vzgadata = data(:,612)/kpc;

% plot galaxy paths
plot3(xgdata, ygdata, zgdata, 'b-');
hold on
plot3(xgadata, ygadata, zgadata, 'g-'); 
axis([-600 600 -800 800 -500 500]);

%animation station
p1 = plot3(xdata(1,:), ydata(1,:), zdata(1,:), 'r*'); % milky way stars
p2 = quiver3(xdata(1,:), ydata(1,:), zdata(1,:), ...
   vxdata(1,:), vydata(1,:), vzdata(1,:), 'b-');
p3 = plot3(xgdata(1), ygdata(1), zgdata(1), 'k+'); %milky way
p4 = quiver3(xgdata(1), ygdata(1), zgdata(1), ...
    vxgdata(1), vygdata(1), vzgdata(1), 'c-');
p5 = plot3(xadata(1,:), yadata(1,:), zadata(1,:), 'g*'); %andromeda stars
p6 = quiver3(xadata(1,:), yadata(1,:), zadata(1,:), ...
    vxadata(1,:), vyadata(1,:), vzadata(1,:), 'm-');
p7 = plot3(xgadata(1), ygadata(1), zgadata(1), 'k+'); %andromeda
p8 = quiver3(xgadata(1), ygadata(1), zgadata(1), ...
    vxgadata(1), vygadata(1), vzgadata(1), 'b-');
plot3(Bx, By, Bz, 'k+', 'MarkerSize', 15); % barycenter

for i=2:size(xdata(:,1))
    set(p1, 'XData', xdata(i,:), 'YData', ydata(i,:), 'ZData', zdata(i,:));
    set(p2, 'XData', xdata(i,:), 'YData', ydata(i,:), ...
        'ZData', zdata(i,:), 'UData', vxdata(i,:), ...
        'VData', vydata(i,:), 'WData', vzdata(i,:));
    set(p3, 'XData', xgdata(i), 'YData', ygdata(i), 'ZData', zgdata(i));
    set(p4, 'XData', xgdata(i), 'YData', ygdata(i), 'ZData', zgdata(i), ...
        'UData', vxgdata(i), 'VData', vygdata(i), 'WData', vzgdata(i));
    set(p5, 'XData', xadata(i,:), 'YData', yadata(i,:), ...
        'ZData', zadata(i,:));
    set(p6, 'XData', xadata(i,:), 'YData', yadata(i,:), ...
        'ZData', zadata(i,:), 'UData', vxadata(i,:), ...
        'VData', vyadata(i,:), 'WData', vzadata(i,:));
    set(p7, 'XData', xgadata(i), 'YData', ygadata(i), 'ZData', zgadata(i));
    set(p8, 'XData', xgadata(i), 'YData', ygadata(i), ...
        'ZData', zgadata(i), 'UData', vxgadata(i), ...
        'VData', vygadata(i), 'WData', vzgadata(i));
    pause(0.001);
end

%% functions 
function [d] = dist(x1, y1, z1, x2, y2, z2)
d = sqrt((x2 - x1).^2 + (y2 - y1).^2 + (z2 - z1).^2);
end
