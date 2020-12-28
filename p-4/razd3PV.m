Dp = 8*10^-3;
Fpu = 20 / (10^-6);
Vp = 170/60;
a = 1000*10^-3;
n0 = 1000;
Dn = 0.83;
Dv = 0.45;
Tokr = 50;
%%
%razd3.2
Nvitkiv=floor(a/Dp);
Nshariv=floor((Dn-Dv)/(2*Dp));
%
l1shar=pi*(Dv+Dp);
l23shar=pi*(Dv+Dp*(2*Nshariv-1));

t1shar=l1shar/Vp*Nvitkiv;
t23shar=l23shar/Vp*Nvitkiv;
tishar=pi/Vp*Nvitkiv;
