clc;
clear;
close all;


global n0 G v

n0=.1; G=64*10^(-6); v=20;

dt=.1;t_t=10;t=0:dt:t_t;
fintered_noise=q_filter(dt,t_t);

plot(t,fintered_noise*100)
ylabel('Road roughness process (cm)')
xlabel('Time (s)')

% success.