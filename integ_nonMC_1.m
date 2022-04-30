clc;
clear;
close all;


% this script gives MCS of **nonlinear** vehicle
% nonlinearity on this script: Duffing tires and bilinear suspensions


global ms mb mp mf mr ks kf kr cs cf cr ktf ktr l1 l2 l3 qf1 qr1 t1  e1 e2 e3 e4 e5

e1=0; e2=0; e3=0; e4=40; e5=40;

ms=100; mb=12000; mp=3500; mf=180; mr=300; 
ks=20000; kf=120000; kr=120000; 
cs=550; cf=1200; cr=1200; 
ktf=520000; ktr=520000; 
l1=.5; l2=1; l3=2;

z0=zeros(10,1);

global n0 G v

n0=.1; G=4*10^(-6); v=15;

dt=.1;t_t=50;t1=0:dt:t_t;

tic
MCn = 500;
for kk = 1:MCn
    qf1=q_filter(dt,t_t); qr1=qf1;
    [t,z]=ode45(@linear_integer,t1,z0);
    ys_mtx(:,kk)=z(:,1);
    yb_mtx(:,kk)=z(:,2);
    ya_mtx(:,kk)=z(:,3);
    yf_mtx(:,kk)=z(:,4);
    yr_mtx(:,kk)=z(:,5);
    
end
toc


for ii = 1:length(t)
    var_s(ii) = var(ys_mtx(ii,:));
    var_b(ii) = var(yb_mtx(ii,:));
    var_a(ii) = var(ya_mtx(ii,:));
    var_f(ii) = var(yf_mtx(ii,:));
    var_r(ii) = var(yr_mtx(ii,:));
end

figure(1)
plot(t,sqrt(var_s))

function dz=linear_integer(t,z)

global ms mb mp mf mr ks kf kr cs cf cr ktf ktr l1 l2 l3 qf1 qr1 t1 e1 e2 e3 e4 e5

qf=interp1(t1,qf1,t);
qr=interp1(t1,qr1,t);

M = [ms,0,0,0,0;0,mb,0,0,0;0,0,mp,0,0;0,0,0,mf,0;0,0,0,0,mr];
C = [cs,-cs,cs*l1,0,0;...
    -cs,cs,-cs*l1,0,0;...
    cs*l1,-cs*l1,cs*l1^2+cr*l3^2,0,0;...
    0,0,0,0,0;...
    0,0,0,0,0];

K = [ks,-ks,ks*l1,0,0;...
    -ks,ks+kf+kr,-ks*l1-kf*l2+kr*l3,-kf,-kr;...
    ks*l1,-ks*l1-kf*l2+kr*l3,ks*l1^2+kf*l2^2+kr*l3^2,kf*l2,-kr*l3;...
    0,-kf,kf*l2,kf+ktf,0;...
    0,-kr,-kr*l3,0,kr+ktr];

Factor_mtx = [0,0;0,0;0,0;ktf,0;0,ktr];

Q = [qf;qr];

non_vec = [0;...
    cf*abs(z(7)-l2*z(8)-z(9));...
    cr*abs(z(7)+l3*z(8)-z(10))*l3-cf*abs(z(7)-z(8)*l2-z(9))*l2;...
    -cf*abs(z(7)-l2*z(8)-z(9))+e4*ktf*(z(4)-qf)^3;...
    -cf*abs(z(7)+l3*z(8)-z(10))+e5*ktr*(z(5)-qr)^3];


vel_vec = [z(6);z(7);z(8);z(9);z(10)];disp_vec = [z(1);z(2);z(3);z(4);z(5)];


dz = [vel_vec;inv(M)*(-non_vec-C*vel_vec-K*disp_vec+Factor_mtx*Q)];




end