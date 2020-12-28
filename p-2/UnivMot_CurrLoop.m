%DP-12, D-12
% control object K / (Ta * s + 1)
% FO - mu = 1.5
clear
Pn=2500;
nn=1100;
Un=220;
In=16;
Rj=1.63;
N=990;
a=1;
p=2;
Rv=0.59;
wv=83;
Fn=0.0052;
nmax=3600;
Jd=0.05;

%Ff=[0.0001,0.29,0.51,0.7,0.82,0.915,0.99,1.05,1.1,1.15,1.2,1.24,1.28].*Fn*2*p*wv;
Ff=[0,0.29,0.51,0.7,0.82,0.915,0.99,1.05,1.1,1.15,1.2,1.24,1.28].*Fn;
If=[0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4].*16;

Ff1=[0,0.08,0.152,0.215,0.27,0.315,0.355,0.39,0.42,0.445,0.47,0.49,0.51,0.525,0.54,0.551,0.56,0.57,0.579,0.587,0.595,0.603,0.61,0.617,0.623,0.627,0.632,0.637,0.64,0.644,0.647, 0.65 ,0.666].*Fn;
If1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 4.0].*In;

kk=p*N/2/pi/a;
Mn=kk*Fn*In;
omega_n=2*pi*nn/60;
La=0.25*Un/In/Rj*60/2/p/pi/nn;
Rsum=(Rj+Rv)*(1+0.004*(180-20));
Ta=La/Rsum;
Tv=Fn*2*p*wv/16;
Jsum=4*Jd;
Tm=Jsum*Rsum/(kk*Fn);

kp=Un/10;
Tmu=0.001;
Kdt=10/(In*3);
Kds=10/(2*pi*nmax/60);
K = 1 / Rsum;

Tsum = Ta+Tv;

dt = 1e-4;
n_point = 8000;

mu_c = 1.5;

k_riem1(n_point) = 0;
x1=abs(0.5);
for j=1:n_point
k_riem1(j)=(dt^x1)/gamma(x1)*(j^x1-(j-1)^x1)/x1;
end

a = 0.37;
b=exp(-10.27+7.831*mu_c)/a;
% a = 0.4;
% b = 6;
K1 = Tsum^(1-mu_c) * kp * Kdt * K / a;
K2 = Tsum^(-mu_c) * kp * Kdt * K / (a * b);


ua(n_point) = 0;
de(n_point) = 0;
dpsi(n_point) = 0;
ia(n_point) = 0;
fi(n_point) = 0;
psi(n_point) = 0;

u_zt(n_point) = 0;
du_t(n_point) = 0;
u_1(n_point) = 0;
u_21(n_point) = 0;
u_22(n_point) = 0;
u_2(n_point) = 0;
u_sumR(n_point) = 0;
u_dt(n_point) = 0;

for i=1:n_point
    % reference current
    u_zt(i) = 8;
    
    % summator
    if i==1
        du_t(i) = u_zt(i)-0;
    else
        du_t(i) = u_zt(i)-u_dt(i-1);
    end
    
    %-1- I^0.5 -----------------------------------
    s=0;
    for j=1:i
        s = s + du_t(i-j+1)*k_riem1(j);
    end

    u_1(i)= K1 * s;
    
    %-2- I^2 -----------------------------------
    if i==1
        u_21(i)= 0 + du_t(i)*1*dt;
    else
        u_21(i)= u_21(i-1) + du_t(i)*1*dt;
    end
    
    if i==1
        u_22(i)= 0 + u_21(i)*1*dt;
    else
        u_22(i)= u_22(i-1) + u_21(i)*1*dt;
    end
    
    u_2(i) = K2 * u_22(i);
    
    % reg sum
    u_sumR(i) = u_1(i) + u_2(i);
    
%     if u_sumR(i) > 10 %#
%         u_sumR(i) = 10;
%     elseif u_sumR(i) < -10
%         u_sumR(i) = -10;
%     end
    
    % power converter
    ua(i) = kp * u_sumR(i);
    
    if i==1
        dpsi(i) = ua(i)-0;
    else
        dpsi(i) = ua(i)-psi(i-1);
    end
    
      % nonlinear
%     if i==1
%         ia(i)= 0.5*((0*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (0*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
%     else
%         ia(i)= 0.5*((ia(i-1)*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (ia(i-1)*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
%     end
    
    % linear
    if i==1
        ia(i)= 0.5*((0*Tsum + 1 / Rsum * ua(i)*dt)/ (Tsum + dt) + (0*(Tsum - dt) + 1 / Rsum * ua(i)*dt)/Tsum);
    else
        ia(i)= 0.5*((ia(i-1)*Tsum + 1 / Rsum * ua(i)*dt)/ (Tsum + dt) + (ia(i-1)*(Tsum - dt) + 1 / Rsum * ua(i)*dt)/Tsum);
    end
    ia(i) = ia(i) + (randi(10)-5)*0.015;
    
    % nonlinear function of FLUX
    fi(i) = spline(If1,Ff1,ia(i));
    
    % fluxlinkage
    if i==1
         psi(i) = 4*p*wv*(fi(i) - 0) / dt;
    else
         psi(i) = 4*p*wv*(fi(i) - fi(i-1)) / dt;
    end
    
    % feedback current
    u_dt(i)= Kdt *ia(i);
    
     
end
