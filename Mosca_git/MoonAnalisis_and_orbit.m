clear all
clc 
format long

% Kernels loading
cspice_furnsh( '.\kernels\mykernels.furnsh' );

rur=cspice_bodvrd( 'Uranus', 'RADII', 3 );
rur=rur(2);%[km]

gm_ur=cspice_bodvrd('Uranus','GM',1);
gm_sun=cspice_bodvrd('Sun','GM',1);

d_ur_sun=2870660000;%[km]

rSoI=600000;%[km]
rSoI=[rSoI;rSoI];


date0 = '2008 OCT 30 17:14:37.68'; %t0 in UTC
ET_START = cspice_str2et(date0); %t0 in ET

[x_titania_0,ltime]= cspice_spkezr('TITANIA',ET_START,'J2000','NONE','URANUS');

[x_oberon_0,ltime]= cspice_spkezr('OBERON',ET_START,'J2000','NONE','URANUS');

%integrazione

etc1=[ET_START:10:ET_START+1277600];

options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-14);

[time,state_titania] = ode113(@planetary_trajectory,etc1,x_titania_0,options);
[time,state_oberon] = ode113(@planetary_trajectory,etc1,x_oberon_0,options);

%plot
whitebg('black')

plot3(state_titania(:,1),state_titania(:,2),state_titania(:,3),'r')
axis equal
axis([-rSoI(2,1) rSoI(2,1) -rSoI(2,1) rSoI(2,1) -rSoI(2,1)  rSoI(2,1)])
hold on

plot3(state_oberon(:,1),state_oberon(:,2),state_oberon(:,3),'b')
axis equal
axis([-rSoI(2,1) rSoI(2,1) -rSoI(2,1) rSoI(2,1) -rSoI(2,1)  rSoI(2,1)])
hold on

grid on

% PLOT Urano

[x_plane_eq,y_plane_eq]=meshgrid(-rSoI:70000:rSoI);

plot3(x_plane_eq,y_plane_eq,zeros(length(x_plane_eq),length(x_plane_eq)),'w')

[xs,ys,zs] = sphere;
xs = rur(1,1)*xs;
ys = rur(1,1)*ys;   
zs = rur(1,1)*zs;
surf(xs,ys,zs)

sph1 = findobj('Type', 'surface');
cspice_kclear % Plot Uranus surface on the sphere
im1 = imread('uranus.jpg');
set(sph1, 'CData', im1, 'FaceColor', 'texturemap', 'FaceAlpha', 1);

grid on



