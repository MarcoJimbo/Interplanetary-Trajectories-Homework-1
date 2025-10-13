function [rsp_body,x_body,grav_body] = gravity_3(body,frame,system_center,sp_pos,t)

[x_body,ltime]= cspice_spkpos(body,t,frame,'NONE',system_center);

rsp_body = sqrt((sp_pos(1)-x_body(1))^2+(sp_pos(2)-x_body(2))^2+(sp_pos(3)-x_body(3))^2);

GM = cspice_bodvrd( body, 'GM', 3);
%GM= 8.978138845307376e+03;


grav_body = -GM/rsp_body^3;