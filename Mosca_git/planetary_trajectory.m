function dx = planetary_trajectory(t,x)

%##########################################################################
% I. Designate Global Structure
%##########################################################################
  
%##########################################################################
% II. Dynamical Equations d(X) = A X + B
%##########################################################################

%##########################################################################
% II.1 Gravity initialization
% The function gravity determines the position of the 
% spacecraft relative to a generic body(-GM/rsp^3)
%##########################################################################

[rsp,pcb,grav] = gravity_3('Uranus','J2000','Uranus',x(1:3),t);

ID = cspice_bodn2c('Uranus');

% Matrix A

A =[0,0,0,1,0,0
    0,0,0,0,1,0
    0,0,0,0,0,1
    grav,0,0,0,0,0
    0,grav,0,0,0,0
   0,0,grav,0,0,0];

dx = A*x; 
      
      