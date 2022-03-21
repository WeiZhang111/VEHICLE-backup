clc;
clear;
close all;

G=64*10^(-6);v=16.67;
w=0.1:.01:100;
PSD1=G*v./(w.^2+9);
plot(w,PSD1)

xlim([3,100])