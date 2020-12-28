%DP-12, D-12
% 
% 
clear
Pn=2500;
nn=1100;
Un=220;
In=16;
Rj=1.43;
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

Ff1=[0,0.08,0.152,0.215,0.27,0.315,0.355,0.39,0.42,0.445,0.47,0.49, ...
     0.51,0.525,0.54,0.551,0.56,0.57,0.579,0.587,0.595,0.603,0.61, ...
     0.617,0.623,0.627,0.632,0.637,0.64,0.644,0.647, 0.65 ,0.666].*2.2 *Fn;
If1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, ...
     1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, ...
     2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 4.0].*In;

kk=p*N/2/pi/a;
Mn=kk*Fn*In;
omega_n=2*pi*nn/60;
La=0.25*Un/In/Rj/2/p/pi/nn;
Rsum=(Rj+Rv)*(1+0.004*(180-20));
Ta=La/Rsum;
Tv=Fn*2*p*wv/16;
Jsum=4*Jd;
Tm=Jsum*Rsum/(kk*Fn);

kp=Un/10;

Kdt=10/(In*2);
Kds=10/(2*pi*nn/60);
K = 1 / Rsum;
Tmu=0.115;
Tsum = Ta+Tv;

dt = 1e-5;
n_point = 20000;

mu_c = 1.5;

k_riem1(n_point) = 0;
x1=abs(0.5);
for j=1:n_point
k_riem1(j)=(dt^x1)/gamma(x1)*(j^x1-(j-1)^x1)/x1;
end

% current coeff
% a = 0.33;
% % b=exp(-10.27+7.831*mu_c)/a; a = 0.4;
% b = 8;
% K0_c = 1 / (a * b * Tmu^mu_c * kp * Kdt * K);
% K1_c = K0_c * b * Tmu * Tsum; % D^0.5
% K2_c = K0_c * (b * Tmu + Tsum); % I^0.5
% K3_c = K0_c; % I^1.5

K_c = 3.408;
T1 = 0.0006032;
T0 = 0.04728;
mu = 0.5;
Te = 0.2;

znam = T1*dt + T0*dt*dt;
hnK = K_c*dt*dt/znam;
hnA1 = (T1*dt)/znam;
hnB0 = (1*dt*dt)/znam;

% speed coeff
% a = 0.2;
% b=exp(-10.27+7.831*mu_c)/a;
% b = 10;

% K0_s = Jsum / (a * b * Te^mu_c * kk * Fn * Kds * K_c);
% K1_s = b * K0_s * Te * T1; % D^2
% K2_s = K0_s * (b * Te * T0 + T1); % D
% K3_s = K0_s * b * Te; % D^0.5
% K4_s = K0_s * T0; % P
% K5_s = K0_s; % I^1.5



% ua(n_point) = 0;
% de(n_point) = 0;
% dpsi(n_point) = 0;
% ia(n_point) = 0;
% fi(n_point) = 0;
% psi(n_point) = 0;
% 
% u_zt(n_point) = 0;
% du_t(n_point) = 0;
% u_1c(n_point) = 0;
% u_2c(n_point) = 0;
% u_sumR(n_point) = 0;
% u_dt(n_point) = 0;
dy = 0;
for i=1:n_point
    u_zs(i) = 8;
    
    % speed summator
    if i==1
        du_s(i) = u_zs(i);
    else
        du_s(i) = u_zs(i)-u_ds(i-1);
    end
    
    % speed regulator 
    %-1- D^2 --------------------------------
    if i==1
        u_1s(i) = K1_s*(du_s(i) - 2 * 0 + 0)/dt^2;
    elseif i==2
        u_1s(i) = K1_s*(du_s(i) - 2 * du_s(i-1) + 0)/dt^2;
    else
        u_1s(i) = K1_s*(du_s(i) - 2 * du_s(i-1) + du_s(i-2))/dt^2;
    end
    
%     if u_1s(i) > 10 %#
%         u_1s(i) = 10;
%     elseif u_1s(i) < -10
%         u_1s(i) = -10;
%     end
    
    %-2- D -----------------------------------
    if i==1
        u_2s(i) = K2_s*(du_s(i) - 0)/dt;
    else
        u_2s(i) = K2_s*(du_s(i) - du_s(i-1))/dt;
    end
    
    %-3- D^0.5 --------------------------------
    s=0;
    for j=1:i
        s = s + du_s(i-j+1)*k_riem1(j);
    end
    
    u_33s(i)= s;
    
    if i==1
        u_3s(i) = K3_s*(u_33s(i) - 0)/dt;
    else
        u_3s(i) = K3_s*(u_33s(i) - u_33s(i-1))/dt;
    end
    
    %-4- P --------------------------------
    u_4s(i) = K4_s * du_s(i);
    
    %-5- I^1.5 --------------------------------
    if i==1
        u_5s(i) = K5_s * s * dt;
    else
        u_5s(i) = u_5s(i-1) + K5_s * s * dt;
    end
    

    
    % reg sum
    u_s(i) = u_1s(i) + u_2s(i) + u_3s(i) + u_4s(i) +  u_5s(i);
    
     
    %-0- Filter --------------------------------
    if i==1
         u_0s(i) = 0.5*((0*Te + 1*u_s(i)*dt)/ (Te + dt) + (0*(Te - dt) + 1*u_s(i)*dt)/Te);
    else
         u_0s(i)= 0.5*((u_0s(i-1)*Te + 1*u_s(i)*dt)/ (Te + dt) + (u_0s(i-1)*(Te - dt) + 1*u_s(i)*dt)/Te);
    end
    
%     if u_0s(i) > 10 %#
%         u_0s(i) = 10;
%     elseif u_0s(i) < -10
%         u_0s(i) = -10;
%     end
    
    % reference current
    u_zt(i) = u_0s(i);
    
    
    % ideal    
%     x = u_0s(i);
%     
%     dy0 = dy;
%     dy = hnK*x;
%     if i > 1
%         dy = dy + hnA1*dy0 - hnB0*y;
%     end
%     
%     f1(i) = dy;
%     
%     fmu1(i) = f1(i)*(dt^mu)/gamma(1+mu);
%     Imu1 = 0;
%     for j=1:i
%         Imu1 = Imu1 + fmu1(i-j+1)*(j^mu -(j-1)^mu);
%     end
%     y = Imu1;
%     ia(i) = Imu1;
    
    % nonideal
    % ===============================================================
    
    %current summator
    if i==1
        du_t(i) = u_zt(i)-0;
    else
        du_t(i) = u_zt(i)-u_dt(i-1);
    end
    
    %-1- D^0.5 -----------------------------------
    s=0;
    for j=1:i
        s = s + du_t(i-j+1)*k_riem1(j);
    end
    
    u_12c(i) = s;
    if i==1
        u_1c(i) = K1_c*(u_12c(i) - 0)/dt;
    else
        u_1c(i) = K1_c*(u_12c(i) - u_12c(i-1))/dt;
    end
    
    %-2- I^0.5 -----------------------------------
    u_2c(i)= K2_c * s;
    
    %-3- I^1.5 -----------------------------------    
    if i==1
        u_3c(i) = K3_c * s * dt;
    else
        u_3c(i) = u_3c(i-1) + K3_c * s * dt;
    end
    
    % reg sum
    u_sumR(i) = u_1c(i) + u_2c(i) + u_3c(i);
    
%     if u_sumR(i) > 10 %#
%         u_sumR(i) = 10;
%     elseif u_sumR(i) < -10
%         u_sumR(i) = -10;
%     end

    % power converter    
    if i==1
        ua(i) = 0.5*((0*Tmu + kp*u_sumR(i)*dt)/ (Tmu + dt) + (0*(Tmu - dt) + kp*u_sumR(i)*dt)/Tmu);
    else
        ua(i)= 0.5*((ua(i-1)*Tmu + kp*u_sumR(i)*dt)/ (Tmu + dt) + (ua(i-1)*(Tmu - dt) + kp*u_sumR(i)*dt)/Tmu);
    end
    
    if ua(i) > 220
        ua(i) = 220;
    elseif ua(i) < -220
        ua(i) = -220;
    end
    
    if i==1
        de(i) = 0;
    else
        de(i) = ua(i)-e(i-1);
    end
    
    if i==1
        dpsi(i) = 0;
    else
        dpsi(i) = de(i)-psi(i-1);
    end
    
      % nonlinear
    if i==1
        ia(i)= 0.5*((0*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (0*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
    else
        ia(i)= 0.5*((ia(i-1)*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (ia(i-1)*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
    end
    
    %linear
%     if i==1
%         ia(i) = 0;
%     else
%         ia(i)= 0.5*((ia(i-1)*Tsum + 1 / Rsum * ua(i)*dt)/ (Tsum + dt) + (ia(i-1)*(Tsum - dt) + 1 / Rsum * ua(i)*dt)/Tsum);
%     end
%     ia(i) = ia(i) + (randi(10)-5)*0.015;
    
    % ===============================================================
    
    % nonlinear function of FLUX
    fi(i) = spline(If1,Ff1,abs(ia(i))); 
    
    % fluxlinkage
    if i==1
         psi(i) = 2*p*wv*fi(i) / dt;
    else
         psi(i) = 2*p*wv*(fi(i) - fi(i-1)) / dt;
    end
    
    % feedback current
    u_dt(i)= Kdt *ia(i);
    
    % kFi
    kF = fi(i) * kk;
%     kF = Fn * kk; % constant field
    
    % torque
    M(i) = ia(i) * kF;
    
    % torque load
    if i==1
         Mv = 0;
    else
         Mv = 0.0019759 * w(i-1)^2;
    end
    
    % dinamic torque
     Mdin = M(i) - 0;
%        Mdin = M(i) - 0.2 * Mn;
    
    % speed
    if i==1
        w(i)= Mdin / Jsum * dt;
    else
        w(i)= w(i-1) + Mdin / Jsum * dt;
    end
    
    % back electromotive force
    e(i) = w(i) * kF;
    
    % feedback speed
    u_ds(i)= Kds * w(i);
    
end

(max(w)-u_zs(1)/Kds)/(u_zs(1)/Kds)*100
