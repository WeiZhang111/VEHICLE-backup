clc;
clear;
close all;

global ms mb mp mf mr ks kf kr cs cf cr ktf ktr l1 l2 l3 qf1 qr1 t1

ms=70; mb=2100; mp=3500; mf=140; mr=210; ks=12200; kf=74000; kr=120000; cs=550; cf=180; cr=1200; ktf=520000; ktr=520000; l1=.5; l2=2; l3=2;
z0=zeros(10,1);

% t1=0:.1:10;
% qf1=(1/500)*randn(length(t1),1);qr1=(1/500)*randn(length(t1),1);


global n0 G v

n0=.1; G=64*10^(-6); v=16.6667;

dt=.01;t_t=10;t1=0:dt:t_t;
qf1=q_filter(dt,t_t); qr1=q_filter(dt,t_t);


[t,z]=ode45(@linear_integer,t1,z0);

figure(1)
plot(t,z(:,1))
title('Vertical displacement of the seat (m)')
figure(2)
plot(t,z(:,2))
title('Vertical displacement of vehicle body (m)')
figure(3)
plot(t,z(:,3))
title('Angular displacement of vehicle body (rad)')
figure(4)
plot(t,z(:,4))
title('Vertical displacement of front suspension (m)')
figure(5)
plot(t,z(:,5))
title('Vertical displacement of rear suspension (m)')
figure(6)
plot(t1,qf1)
title('Road roughness (m)')

function dz=linear_integer(t,z)

global ms mb mp mf mr ks kf kr cs cf cr ktf ktr l1 l2 l3 qf1 qr1 t1

qf=interp1(t1,qf1,t);
qr=interp1(t1,qr1,t);

dz(1)=z(6);
dz(2)=z(7);
dz(3)=z(8);
dz(4)=z(9);
dz(5)=z(10);
dz(6)=-(1/ms)*(cs*z(6)-cs*z(7)+cs*l1*z(8)+ks*z(1)-ks*z(2)+ks*l1*z(3));
dz(7)=-(1/mb)*(-cs*z(6)+(cs+cf+cr)*z(7)+(-cs*l1-cf*l2+cr*l3)*z(8)+(-cf)*z(9)+(-cr)*z(10)+(-ks)*z(1)+(ks+kf+kr)*z(2)+(-ks*l1-kf*l2+kr*l3)*z(3)+(-kf)*z(4)+(-kr)*z(5));
dz(8)=-(1/mp)*((cs*l1)*z(6)+(-cs*l1-cf*l2+cr*l3)*z(7)+(cs*l1^2+cf*l2^2+cr*l3^2)*z(8)+(cf*l2)*z(9)+(-cr*l3)*z(10)+(ks*l1)*z(1)+(-ks*l1-kf*l2+kr*l3)*z(2)+(ks*l1^2+kf*l2^2+kr*l3^2)*z(3)+(kf*l2)*z(4)+(-kr*l3)*z(5));
dz(9)=(1/mf)*ktf*qf-(1/mf)*(-cf*z(7)+cf*l2*z(8)+cf*z(9)-kf*z(2)+kf*l2*z(3)+(kf+ktf)*z(4));
dz(10)=(1/mr)*ktr*qr-(1/mr)*(-cr*z(7)-cr*l3*z(8)+cr*z(10)-kr*z(2)-kr*l3*z(3)+(kr+ktr)*z(5));
dz=[dz(1);dz(2);dz(3);dz(4);dz(5);dz(6);dz(7);dz(8);dz(9);dz(10)];
end