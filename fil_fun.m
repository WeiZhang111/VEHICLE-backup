function dq = fil_fun(t,q)

global W_1 tspan_1 n0 G v
Wt=interp1(tspan_1,W_1,t);
dq = 2*pi*.1*sqrt(G*v)*Wt-3*v*q;   % linear


end