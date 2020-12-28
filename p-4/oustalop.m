%function G=new_fod(r,N,wb,wh,b,d)
%if nargin==4, b=10; d=9; end
clear
gam = -0.5;
N = 16;
wb = 5e-4;
wh = 5e6;

Ts = 1e-5;

k=1:N; 
wu=sqrt(wh/wb);
K = wh^gam;
wkp=wb*wu.^((2*k-1+gam)/N); 
wk=wb*wu.^((2*k-1-gam)/N);
Gs=zpk(-wk,-wkp, K);
Gd=tf(Gs);

Hd = c2d(Gd,Ts,'tustin');
zN = Hd.Numerator{1, 1};
zD = Hd.Denominator{1, 1};

% zNr = round(zN, 8);
% zDr = round(zD, 8);

%
% gam = -0.3;
% N = 6;
% wb = 1e-4;
% wh = 1e4;
% b = 10;
% d = 9;
% 
% k=1:N; wu=sqrt(wh/wb);
% K=(d*wh/b)^gam;
% wkp=wb*wu.^((2*k-1+gam)/N); 
% wk=wb*wu.^((2*k-1-gam)/N);
% G=zpk(-wk,-wkp, K)*tf([d,b*wh,0],[d*(1-gam),b*wh,d*gam]);
% G=tf(G);
% 
% Hd = c2d(G,0.01,'tustin')
% zN = Hd.Numerator{1, 1};
% zD = Hd.Denominator{1, 1};