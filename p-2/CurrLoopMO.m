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
% Ff=[0,0.29,0.51,0.7,0.82,0.915,0.99,1.05,1.1,1.15,1.2,1.24,1.28].*Fn;
% If=[0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4].*16;

Ff1=[0,0.08,0.152,0.215,0.27,0.315,0.355,0.39,0.42,0.445,0.47,0.49, ...
     0.51,0.525,0.54,0.551,0.56,0.57,0.579,0.587,0.595,0.603,0.61, ...
     0.617,0.623,0.627,0.632,0.637,0.64,0.644,0.647, 0.65 ,0.666].*2.2*Fn;
If1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, ...
     1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, ...
     2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 4.0].*In;
 If2 = 0:0.01:70;
  a =    0.006728;
  b =    0.002175;
  c =   -0.006706;
  d =    -0.08989;
 Ff2 = a*exp(b*If2) + c*exp(d*If2);
 
 
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
Tmu=0.01;
Tsum = 0.036;

dt = 1e-5;
n_point = 20000;
K1 = Tsum / (2 * Tmu * kp * Kdt * K);
K2 = 1 / (2 * Tmu * kp * Kdt * K);

A = K * dt / (2 * Tsum + dt);
B = -(2 * Tsum - dt) / (2 * Tsum + dt);

C = kp * dt / (2 * Tmu + dt);
D = -(2 * Tmu - dt) / (2 * Tmu + dt);

Tex = 0.0257;

ua(n_point) = 0;
de(n_point) = 0;
dpsi(n_point) = 0;
ia(n_point) = 0;
fi(n_point) = 0;
psi(n_point) = 0;

u_zt(n_point) = 0;
du_t(n_point) = 0;
u_1(n_point) = 0;
u_12(n_point) = 0;
u_3(n_point) = 0;
u_41(n_point) = 0;
u_4(n_point) = 0;

u_21(n_point) = 0;
u_22(n_point) = 0;
u_2(n_point) = 0;
u_sumR(n_point) = 0;
u_dt(n_point) = 0;

max_t = 0.0;
flag = 0;

time_in = 0.0;

for i=1:n_point
    % reference current
    u_zt(i) = 9;
    
    % summator
    if i==1
        du_t(i) = u_zt(i)-0;
    else
        du_t(i) = u_zt(i)-u_dt(i-1);
    end
    
    %-1- P     -----------------------------------    
    u_1(i) = K1 * du_t(i);
    
    %-2- I     -----------------------------------
    if i==1
        u_2(i)= K2 * du_t(i)*dt / 2;
    else
        u_2(i)=  (u_2(i-1) +  K2 * du_t(i)*dt / 2 + K2 * du_t(i-1) * dt / 2);
    end
    
%     if u_2(i) > 10
%         u_2(i) = 10;
%     elseif u_2(i) < -10
%         u_2(i) = -10;
%     end
    
    % reg sum
    u_sumR(i) = u_1(i) + u_2(i);
    
    if u_sumR(i) > 10
        u_sumR(i) = 10;
    elseif u_sumR(i) < -10
        u_sumR(i) = -10;
    end
   
    % power converter
    if i==1
        ua(i) = C * u_sumR(i);
    else
        ua(i)= C * u_sumR(i) + C * u_sumR(i-1) - D * ua(i-1);
    end

    % linear
    if i==1
        ia(i)= A * ua(i);
    else
        ia(i)= A * ua(i) + A * ua(i-1) - B * ia(i-1);
    end
    %ia(i) = ia(i) + (randi(10)-5)*0.015;
    
    % feedback current
    u_dt(i)= Kdt *ia(i);
    
    if i > 1 && flag < 1
        if u_dt(i)>u_dt(i-1)
            max_t = (i-1)*dt;
        else
            flag = 1;
        end
    end
    
    if i > 1
        if u_dt(i) < u_zt(1) - u_zt(1)*0.02
            time_in = (i-1)*dt;
        end
    end
    
end

ton = u_zt(1)*0.02;
up_ton = u_zt(1) + ton;
down_ton = u_zt(1) - ton;
maximum = max(u_dt);
up_ton = up_ton / Kdt;
down_ton = down_ton / Kdt;
maximum = maximum / Kdt;
plot([0 (n_point-1)*dt],[up_ton up_ton],[0 (n_point-1)*dt],[down_ton down_ton],[0 max_t],[maximum maximum],[time_in time_in],[0 down_ton],(0:n_point-1)*dt,ia)


(max(ia)-u_zt(1)/Kdt)/(u_zt(1)/Kdt)*100
%%
%NONLINEAR!!
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


Ff1=[0,0.08,0.152,0.215,0.27,0.315,0.355,0.39,0.42,0.445,0.47,0.49, ...
     0.51,0.525,0.54,0.551,0.56,0.57,0.579,0.587,0.595,0.603,0.61, ...
     0.617,0.623,0.627,0.632,0.637,0.64,0.644,0.647, 0.65 ,0.666].*2.2*Fn;
If1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, ...
     1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, ...
     2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 4.0].*In;
 If2 = 0:0.01:70;
  a =    0.006728;
  b =    0.002175;
  c =   -0.006706;
  d =    -0.08989;
 Ff2 = a*exp(b*If2) + c*exp(d*If2);
 
 
pp = polyfit(If1,Ff1,9);
 
kk=p*N/2/pi/a;
Mn=kk*Fn*In;
omega_n=2*pi*nn/60;
La=0.25*Un/In*60/2/p/pi/nn;
Rsum=(Rj+Rv)*(1+0.004*(180-20));
Ta=La/Rsum;
Jsum=4*Jd;
Tm=Jsum*Rsum/(kk*Fn);

kp=Un/10;

Kdt=10/(In*2);
Kds=10/(2*pi*nn/60);
K = 1 / Rsum;
Tmu=0.01;
Tsum = 0.036;

dt = 1e-5;
n_point = 20000;
K1 = Tsum / (2 * Tmu * kp * Kdt * K);
K2 = 1 / (2 * Tmu * kp * Kdt * K);

A = K * dt / (2 * Ta + dt);
B = -(2 * Ta - dt) / (2 * Ta + dt);

C = kp * dt / (2 * Tmu + dt);
D = -(2 * Tmu - dt) / (2 * Tmu + dt);

Tex = 0.0257;

ua(n_point) = 0;
de(n_point) = 0;
dpsi(n_point) = 0;
ia(n_point) = 0;
fi(n_point) = 0;
psi(n_point) = 0;

u_zt(n_point) = 0;
du_t(n_point) = 0;
u_1(n_point) = 0;
u_12(n_point) = 0;
u_3(n_point) = 0;
u_41(n_point) = 0;
u_4(n_point) = 0;

u_21(n_point) = 0;
u_22(n_point) = 0;
u_2(n_point) = 0;
u_sumR(n_point) = 0;
u_dt(n_point) = 0;

max_t = 0.0;
flag = 0;

time_in = 0.0;

for i=1:n_point
    % reference current
    u_zt(i) = 5;
    
    % summator
    if i==1
        du_t(i) = u_zt(i)-0;
    else
        du_t(i) = u_zt(i)-u_dt(i-1);
    end
    
    %-1- P     -----------------------------------    
    u_1(i) = K1 * du_t(i);
    
    %-2- I     -----------------------------------
    if i==1
        u_2(i)= K2 * du_t(i)*dt / 2;
    else
        u_2(i)=  (u_2(i-1) +  K2 * du_t(i)*dt / 2 + K2 * du_t(i-1) * dt / 2);
    end
    
%     if u_2(i) > 10
%         u_2(i) = 10;
%     elseif u_2(i) < -10
%         u_2(i) = -10;
%     end
    
    % reg sum
    u_sumR(i) = u_1(i) + u_2(i);
    
    if u_sumR(i) > 10
        u_sumR(i) = 10;
    elseif u_sumR(i) < -10
        u_sumR(i) = -10;
    end
   
    % power converter
    if i==1
        ua(i) = C * u_sumR(i);
    else
        ua(i)= C * u_sumR(i) + C * u_sumR(i-1) - D * ua(i-1);
    end
    
    ua(i) = 90;
    if i==1
        dpsi(i) = ua(i)-0;
    else
        dpsi(i) = ua(i)- psi(i-1);
    end
    
      % nonlinear
    if i==1
        ia(i)= A * dpsi(i);
    else
        ia(i)= A * dpsi(i) + A * dpsi(i-1) - B * ia(i-1);
    end

    %ia(i) = ia(i) + (randi(10)-5)*0.015;
    
    % nonlinear function of FLUX
%     if(ia(i) > 0.1)
%         fi(i) = interp1(If1,Ff1,ia(i), 'spline');
%     else
%         fi(i) = 0;
%     end
    
    fi(i) = abs(interp1(If2,Ff2,abs(ia(i)), 'linear'));

%     fi(i) = polyval(pp,ia(i));
%     fi(i) = i*0.00001;
    % fluxlinkage
%     if i==1
%          psi(i) = 2*p*wv*(2 / dt * fi(i));
%     else
%          psi(i) = 2*p*wv*(2 / dt * fi(i) - 2 / dt *fi(i-1) - psi(i-1));
%     end

    if i==1
         psi(i) = 2*p*wv*(fi(i))/dt;
    else
         psi(i) = 2*p*wv*(fi(i) - fi(i-1))/dt;
    end
    
    % feedback current
    u_dt(i)= Kdt *ia(i);
    
    if i > 1 && flag < 1
        if u_dt(i)>u_dt(i-1)
            max_t = (i-1)*dt;
        else
            flag = 1;
        end
    end
    
    if i > 1
        if u_dt(i) < u_zt(1) - u_zt(1)*0.02
            time_in = (i-1)*dt;
        end
    end
    
end

ton = u_zt(1)*0.02;
up_ton = u_zt(1) + ton;
down_ton = u_zt(1) - ton;
maximum = max(u_dt);
up_ton = up_ton / Kdt;
down_ton = down_ton / Kdt;
maximum = maximum / Kdt;
plot([0 (n_point-1)*dt],[up_ton up_ton],[0 (n_point-1)*dt],[down_ton down_ton],[0 max_t],[maximum maximum],[time_in time_in],[0 down_ton],(0:n_point-1)*dt,ia)


(max(ia)-u_zt(1)/Kdt)/(u_zt(1)/Kdt)*100