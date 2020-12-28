%DP-12, D-12
% control object kp / (Tmu * s + 1)
% FO - mu = 1.5
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

Ff1=[0,0.08,0.152,0.215,0.27,0.315,0.355,0.39,0.42,0.445,0.47,0.49,0.51,0.525,0.54,0.551,0.56,0.57,0.579,0.587,0.595,0.603,0.61,0.617,0.623,0.627,0.632,0.637,0.64,0.644,0.647, 0.65 ,0.666].*Fn;
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
Tmu=0.01;
Tsum = 0.036;

dt = 1e-5;
n_point = 20000;

mu_c = 1.5;

k_riem1(n_point) = 0;
x1=abs(0.5);
for j=1:n_point
k_riem1(j)=(dt^x1)/gamma(x1)*(j^x1-(j-1)^x1)/x1;
end

k_riem2(n_point) = 0;
x1=abs(1.5);
for j=1:n_point
k_riem2(j)=(dt^x1)/gamma(x1)*(j^x1-(j-1)^x1)/x1;
end

a = 0.3;
b=exp(-10.27+7.831*mu_c)/a;
% a = 0.4;
% b = 20;
K0 = 1 / (a * b * Tmu^mu_c * kp * Kdt * K);
K1 = K0 * b * Tmu * Tsum; % D^0.5
K2 = K0 * (b * Tmu + Tsum); % I^0.5
K3 = K0; % I^1.5

% FO I^0.5 approx
%------------------------------
gam = -0.5;
N = 16;
wb = 5e-4;
wh = 5e6;

Ts = dt;

k=1:N; 
wu=sqrt(wh/wb);
F = wh^gam;
wkp=wb*wu.^((2*k-1+gam)/N); 
wk=wb*wu.^((2*k-1-gam)/N);
%------------------------------

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
    
    %-1- D^0.5 -----------------------------------
    s=0;
    for j=1:i
        s = s + du_t(i-j+1)*k_riem1(j);
    end
    
    u_12(i) = s;
    if i==1
        u_1(i) = K1*(u_12(i) - 0)/dt;
    else
        u_1(i) = K1*(u_12(i) - u_12(i-1))/dt;
    end
    
    %-2- I^0.5 -----------------------------------
    u_2(i)= K2 * s;
    
    %-3- I^1.5 -----------------------------------    
    if i==1
        u_3(i) = K3 * s * dt;
    else
        u_3(i) = u_3(i-1) + K3 * s * dt;
    end
        

    % reg sum
    u_sumR(i) = u_1(i) + u_2(i) + u_3(i);
           
    % equivalent
%     if i==1
%         ieq(i)= 0;
%     else
%         ieq(i)= (ieq(i-1)*Tmu*0.5 + 1 / Kdt * u_zt(i)*dt)/(Tmu*0.5+dt);
%     end

    
    % power converter
    if i==1
        ua(i) = 0.5*((0*Tmu + kp*u_sumR(i)*dt)/ (Tmu + dt) + (0*(Tmu - dt) + kp*u_sumR(i)*dt)/Tmu);
    else
        ua(i)= 0.5*((ua(i-1)*Tmu + kp*u_sumR(i)*dt)/ (Tmu + dt) + (ua(i-1)*(Tmu - dt) + kp*u_sumR(i)*dt)/Tmu);
    end
    
%     if ua(i) > 220
%         ua(i) = 220;
%     elseif ua(i) < -220
%         ua(i) = -220;
%     end
   
%     if i==1
%         dpsi(i) = ua(i);
%     else
%         dpsi(i) = ua(i)-psi(i-1);
%     end
    
      % nonlinear
%     if i==1
%         ia(i)= 0.5*((0*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (0*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
%     else
%         ia(i)= 0.5*((ia(i-1)*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (ia(i-1)*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
%     end
    
    % linear
    if i==1
        ia(i)= 0;
    else
        ia(i)= 0.5*((ia(i-1)*Tsum + 1 / Rsum * ua(i)*dt)/ (Tsum + dt) + (ia(i-1)*(Tsum - dt) + 1 / Rsum * ua(i)*dt)/Tsum);
    end
%     ia(i) = ia(i) + (randi(10)-5)*0.015;
    
    % nonlinear function of FLUX
%     fi(i) = spline(If1,Ff1,ia(i));
    
    % fluxlinkage
%     if i==1
%          psi(i) = 4*p*wv*(fi(i) - 0) / dt;
%     else
%          psi(i) = 4*p*wv*(fi(i) - fi(i-1)) / dt;
%     end
    
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
