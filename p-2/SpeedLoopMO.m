%DP-12, D-12
% control object kp / (Tmu * s + 1)
% modular optimum
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

Ff1=[0,0.08,0.152,0.215,0.27,0.315,0.355,0.39,0.42,0.445,0.47,0.49,0.51,0.525,0.54,0.551,0.56,0.57,0.579,0.587,0.595,0.603,0.61,0.617,0.623,0.627,0.632,0.637,0.64,0.644,0.647, 0.65 ,0.666].*2.2*Fn;
If1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 4.0].*In;

kk=p*N/2/pi/a;
Mn=kk*Fn*In;
omega_n=2*pi*nn/60;
La=0.25*Un/In*60/2/p/pi/nn;
Rsum=(Rj+Rv)*(1+0.004*(180-20));
Ta=La/Rsum;
Tv=Fn*2*p*wv/16;
Jsum=4*Jd;
Tm=Jsum*Rsum/(kk*Fn);

kp=Un/10;

Kdt=10/(In*2);
Kds=10/(2*pi*nn/60);
K = 1 / Rsum;
Tmu=0.03;
Tsum = Ta+Tv;

dt = 5e-4;
n_point = 8000;

% current coeff
K1_c = Tsum / (2 * Tmu * kp * Kdt * K); % P
K2_c = 1 / (2 * Tmu * kp * Kdt * K); % I

% speed coeff
K1_s = Jsum * Kdt / ( 2 * (2 * Tmu) * kk * Fn * Kds); % P
K2_s = Jsum * Kdt / ( 8 * (2 * Tmu)^2 * kk * Fn * Kds); % I

% K1_s = 4.522250;
% K2_s = 0.411148;

% K1_s = 2.56443022*1.5;
% K2_s = 0.8821*2.5;

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
w(n_point) = 0;

for i=1:n_point
    u_zs(i) = 8;
    
%     if i > 2666 && i <= 5332
%         u_zs(i) = 4;
%     elseif i > 5332 && i <= 8000
%         u_zs(i) = 8;
%     end
    
    % speed summator
    if i==1
        du_s(i) = u_zs(i);
    else
        du_s(i) = u_zs(i)-u_ds(i-1);
    end
    
    % speed regulator
    %-1- P     -----------------------------------    
    u_1s(i) = K1_s * du_s(i);
    
    
    %-2- I     -----------------------------------
    if i==1
        u_2s(i)= du_s(i)*K2_s*dt;
    else
        u_2s(i)= u_2s(i-1) + du_s(i)*K2_s*dt;
    end
    
    if u_2s(i) > 10 %#
        u_2s(i) = 10;
    elseif u_2s(i) < -10
        u_2s(i) = -10;
    end
    
    u_s(i) = u_1s(i) + u_2s(i);
    
    if u_s(i) > 10 %#
        u_s(i) = 10;
    elseif u_s(i) < 0
        u_s(i) = 0;
    end
    
    % reference current
    u_zt(i) = u_s(i);
    
    
    % ideal
%     if i==1
%         ia(i)= 0.5*((0*2*Tmu + 1 / Kdt * u_zt(i)*dt)/ (2*Tmu + dt) + (0*(2*Tmu - dt) + 1 / Kdt * u_zt(i)*dt)/(2*Tmu));
%     else
%         ia(i)= 0.5*((ia(i-1)*2*Tmu + 1 / Kdt * u_zt(i)*dt)/ (2*Tmu + dt) + (ia(i-1)*(2*Tmu - dt) + 1 / Kdt * u_zt(i)*dt)/(2*Tmu));
%     end
        
    % nonideal
    % ===============================================================
    
    %current summator
    if i==1
        du_t(i) = u_zt(i)-0;
    else
        du_t(i) = u_zt(i)-u_dt(i-1);
    end
    
    % current regulator
    %-1- P     -----------------------------------    
    u_1c(i) = K1_c * du_t(i);
    
    %-2- I     -----------------------------------
    if i==1
        u_2c(i)= du_t(i)*K2_c*dt;
    else
        u_2c(i)= u_2c(i-1) + du_t(i)*K2_c*dt;
    end
    
    if u_2c(i) > 10
        u_2c(i) = 10;
    elseif u_2c(i) < -10
        u_2c(i) = -10;
    end
    
    % reg sum
    u_sumR(i) = u_1c(i) + u_2c(i);
    
    if u_sumR(i) > 10
        u_sumR(i) = 10;
    elseif u_sumR(i) < -10
        u_sumR(i) = -10;
    end

    % power converter    
    if i==1
        ua(i) = 0.5*((0*Tmu + kp*u_sumR(i)*dt)/ (Tmu + dt) + (0*(Tmu - dt) + kp*u_sumR(i)*dt)/Tmu);
    else
        ua(i)= 0.5*((ua(i-1)*Tmu + kp*u_sumR(i)*dt)/ (Tmu + dt) + (ua(i-1)*(Tmu - dt) + kp*u_sumR(i)*dt)/Tmu);
    end
    
    if i==1
        de(i) = ua(i);
    else
        de(i) = ua(i)-e(i-1);
    end
    
    if i==1
        dpsi(i) = de(i);
    else
        dpsi(i) = de(i)-psi(i-1);
    end
    
      % nonlinear
    if i==1
        ia(i)= 0.5*((0*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (0*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
    else
        ia(i)= 0.5*((ia(i-1)*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (ia(i-1)*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
    end
    
    % linear
%     if i==1
%         ia(i) = 0;
%     else
%         ia(i)= 0.5*((ia(i-1)*Tsum + 1 / Rsum * ua(i)*dt)/ (Tsum + dt) + (ia(i-1)*(Tsum - dt) + 1 / Rsum * ua(i)*dt)/Tsum);
%     end
%     %ia(i) = ia(i) + (randi(10)-5)*0.015;
    
    % ===============================================================
    
    % nonlinear function of FLUX
    fi(i) = spline(If1,Ff1,abs(ia(i)));
    
    % fluxlinkage
    if i==1
         psi(i) = 2*p*wv*(fi(i) - 0) / dt;
    else
         psi(i) = 2*p*wv*(fi(i) - fi(i-1)) / dt;
    end
    
    % feedback current
    u_dt(i)= Kdt *ia(i);
    
    % kFi
    kF = fi(i) * kk;
%     kF = Fn * kk; % constant field
    
    % torque
    M(i) = abs(ia(i)) * kF;
    
    % torque load
    if i==1
         Mv = 0;
    else
         Mv = 0.0019759 * w(i-1)^2;
    end
    
    Mr = 0;
    if i > 2
        if w(i-1) > 0.001
            Mr = 0.8*Mn;
        end
    end
    
    % dinamic torque
%     Mdin = M(i) - Mv;
    Mdin = M(i) - Mv;
    
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
%     u_ds(i) = u_ds(i) + (randi(10)-5)*0.02;
    
end

(max(w)-u_zs(1)/Kds)/(u_zs(1)/Kds)*100