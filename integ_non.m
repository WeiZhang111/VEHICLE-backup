clc;
clear;
close all;


% this script gives solution to **nonlinear** vehicle


global ms mb mp mf mr ks kf kr cs cf cr ktf ktr l1 l2 l3 qf1 qr1 t1  e1 e2 e3 e4 e5

e1=0; e2=-20; e3=-20; e4=1; e5=1;

ms=70; mb=12000; mp=3500; mf=180; mr=300; 
ks=20000; kf=120000; kr=120000; 
cs=550; cf=1200; cr=1200; 
ktf=520000; ktr=520000; 
l1=.5; l2=1; l3=2;

z0=zeros(10,1);

% t1=0:.1:10;
% qf1=(1/500)*randn(length(t1),1);qr1=(1/500)*randn(length(t1),1);



global n0 G v

n0=.1; G=4*10^(-6); v=15;

dt=.01;t_t=20;t1=0:dt:t_t;
qf1=q_filter(dt,t_t); qr1=qf1;


[t,z]=ode45(@linear_integer,t1,z0);

figure(1)
plot(t,z(:,1))
title('Vertical displacement of the seat (m)')
figure(2)
plot(t,z(:,2))
title('Vertical displacement of vehicle body (m)')
% figure(3)
% plot(t,z(:,3))
% title('Angular displacement of vehicle body (rad)')
% figure(4)
% plot(t,z(:,4))
% title('Vertical displacement of front suspension (m)')
% figure(5)
% plot(t,z(:,5))
% title('Vertical displacement of rear suspension (m)')
% figure(6)
% plot(t1,qf1)
% title('Road roughness (m)')

function dz=linear_integer(t,z)

global ms mb mp mf mr ks kf kr cs cf cr ktf ktr l1 l2 l3 qf1 qr1 t1 e1 e2 e3 e4 e5

qf=interp1(t1,qf1,t);
qr=interp1(t1,qr1,t);

M = [ms,0,0,0,0;0,mb,0,0,0;0,0,mp,0,0;0,0,0,mf,0;0,0,0,0,mr];
C = [cs,-cs,cs*l1,0,0;...
    -cs,cs+cf+cr,-cs*l1-cf*l2+cr*l3,-cf,-cr;...
    cs*l1,-cs*l1-cf*l2+cr*l3,cs*l1^2+cf*l2^2+cr*l3^2,cf*l2,-cr*l3;...
    0,-cf,cf*l2,cf,0;...
    0,-cr,-cr*l3,0,cr];

K = [ks,-ks,ks*l1,0,0;...
    -ks,ks+kf+kr,-ks*l1-kf*l2+kr*l3,-kf,-kr;...
    ks*l1,-ks*l1-kf*l2+kr*l3,ks*l1^2+kf*l2^2+kr*l3^2,kf*l2,-kr*l3;...
    0,-kf,kf*l2,kf+ktf,0;...
    0,-kr,-kr*l3,0,kr+ktr];

Factor_mtx = [0,0;0,0;0,0;ktf,0;0,ktr];

Q = [qf;qr];

non_vec = [e1*ks*(z(1)-z(2)+l1*z(3))^3;...
    -e1*ks*(z(1)-z(2)+l1*z(3))^3+e2*kf*(z(2)-z(3)*l2-z(4))^3+e3*kr*(z(2)+z(3)*l3-z(5))^3;...
    e1*ks*(z(1)-z(2)+l1*z(3))^3*l1-e2*kf*(z(2)-z(3)*l2-z(4))^3*l2+e3*kr*(z(2)+z(3)*l3-z(5))^3*l3;...
    -e2*kf*(z(2)-z(3)*l2-z(4))^3+e4*ktf*(z(4)-qf)^3;...
    -e3*kr*(z(2)+z(3)*l3-z(5))^3+e5*ktr*(z(5)-qr)^3];


vel_vec = [z(6);z(7);z(8);z(9);z(10)];disp_vec = [z(1);z(2);z(3);z(4);z(5)];


dz = [vel_vec;inv(M)*(-non_vec-C*vel_vec-K*disp_vec+Factor_mtx*Q)];




end