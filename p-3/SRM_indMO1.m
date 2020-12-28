clear
Rd = 0.68;

% x = [0 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.75 1.8 1.85 1.9 1.95 2 2.05];
% y = [0 76 86 96 140 190 240 300 400 550 1000 1600 3400 7700 10600 13400 16400 19400 27800 38800 49800];


angleON = 7.5;
angleOFF = 22.5;

Jd = 0.0000073 * 2;

Ts = 2e-5;
n_point = 10000;
time = (0:n_point-1)*Ts;

load('D:\cloud\DISSERT\p - 3 - data\matlab\Current=f(L,Q).mat');
load('D:\cloud\DISSERT\p - 3 - data\matlab\Torque=f(i,Q).mat');

Ix(size(Iv,1),size(Iv,2))=0;
Iy(size(Iv,1),size(Iv,2))=0;
for i = 1:size(Iv,1)
    for j = 1:size(Iv,2)
        Ix(i,j) = (j - 1) * 0.5;
        Iy(i,j) = (i - 1) * 0.5;
    end
end
% i = interp2(Ix,Iy,Iv,L,Q)

Mx(size(Mv,1),size(Mv,2))=0;
My(size(Mv,1),size(Mv,2))=0;
for i = 1:size(Mv,1)
    for j = 1:size(Mv,2)
        Mx(i,j) = (j - 1) * 0.5;
        My(i,j) = (i - 1) * 0.5;
    end
end
% M = interp2(Mx,My,Mv,i,Q)

mu = 0.42;

Tmus = 0.005;
kss = 10 / 500;
Kp = 24.2 / 10;
dt = Ts;

k_riem1(n_point) = 0;
x1=abs(mu - 1);
for j=1:n_point
k_riem1(j)=(dt^x1)/gamma(x1)*(j^x1-(j-1)^x1)/x1;
end

w = 0;
alfa = [0 -15 -30 -45];
alfaMod = mod(alfa, 60);
Ia = 0; Ib = 0; Ic = 0; Id = 0;
La = interp2(Ix,Iy,Iv,Ia,alfaMod(1),'spline'); 
Lb = interp2(Ix,Iy,Iv,Ib,alfaMod(2),'spline');
Lc = interp2(Ix,Iy,Iv,Ic,alfaMod(3),'spline');
Ld = interp2(Ix,Iy,Iv,Id,alfaMod(4),'spline');
psi(4) = 0.0;

Mdim = 0;

current(n_point)=0;
speed(n_point)=0;
torque(n_point)=0;

uv(n_point) = 0;
ufrac(n_point) = 0;
u_c = 0;
u_ss = 0;
u01 = 0;
u00 = 0;
u3 = 0;
u2 = 0;
u4 = 0;
u24 = 0;


for i = 1 : n_point
    alfa = alfa + w / pi * 180 * Ts;
   
    alfaMod = mod(alfa, 60);
    
    u_zs = 6;
    
    if i > 20000 && i <= 40000
        u_zs = 6;
    elseif i > 40000 && i <= 60000
        u_zs = 9;
    end
    
    uv(i) = u_zs;
   
    % speed summator
    du_s = u_zs - u_ss;
    
    % change coeff
    a = 2;
%     b = 2;
%     b=exp(-10.27+7.831*mu_c)/a;
    if u_zs >= 8.33 
        b = 2.2;
        
        a1 = 0.00142; %
        a0 = 0.2364516;
        k = 32.46734; %
    elseif u_zs >= 6.67
        b = 2.5;
        
        a1 = 0.001112; %
        a0 = 0.343870;
        k = 41.71000; %
    elseif u_zs >= 5
        b = 3;
        
        a1 = 0.0011305; %
        a0 = 0.404838;
        k = 48.99223; %
    elseif u_zs >= 3.33
        b = 3.55;
        
        a1 = 0.001496; %
        a0 = 0.481290;
        k = 57.7992; %   
    elseif u_zs >= 1.67
        b = 5.2;
        
        a1 = 0.002795; %
        a0 = 0.565483;
        k = 69.2641; %
    elseif u_zs >= 0
        b = 8;
        
        a1 = 0.006661; %
        a0 = 0.50161;
        k = 69.0880; %
    end
    
    K1s = a1 / (a * Tmus * kss * Kp * k);
    K2s = a0 / (a * Tmus * kss * Kp * k);
    K3s = 1 / (a * Tmus * kss * Kp * k);
    
    ufrac(i) = du_s;
    s=0;
    for j=1:i
        s = s + ufrac(i-j+1)*k_riem1(j);
    end
    u00 = u01;
    u01 = s;
    
    u1 = K1s * (u01 - u00) / dt;
    
    u2 = K2s * s;
    
    u3 = u3 + K3s * du_s *dt;
    
    ureg = u1 + u2 + u3;
    
    u_c= 0.5*((u_c * Tmus + Kp * ureg*dt)/ (Tmus + dt) + (u_c*(Tmus - dt) + Kp * ureg*dt)/Tmus);
    
    if u_c > 0
        angleON = 7.5;
        angleOFF = 22.5;
    else
        angleON = 37.5;
        angleOFF = 52.5;
    end
    
    sig = (alfaMod >= angleON) & (alfaMod <= angleOFF);
    Uv = abs(u_c);
%     Uv = 6;
    
    % voltage on phase
    if (sig(1) < 1 && Ia > 0)
        Ua = -24.2;
    elseif (sig(1) > 0)
        Ua = Uv;
    else
        Ua = 0;
    end
    
    if (sig(2) < 1 && Ib > 0)
        Ub = -24.2;
    elseif (sig(2) > 0)
        Ub = Uv;
    else
        Ub = 0;
    end
    
    if (sig(3) < 1 && Ic > 0)
        Uc = -24.2;
    elseif (sig(3) > 0)
        Uc = Uv;
    else
        Uc = 0;
    end
    
    if (sig(4) < 1 && Id > 0)
        Ud = -24.2;
    elseif (sig(4) > 0)
        Ud = Uv;
    else
        Ud = 0;
    end
    
    Uall = [Ua Ub Uc Ud] - [Ia Ib Ic Id] * Rd;
    psi = psi + 1 * Uall * Ts;
    
    Ia = psi(1) / La;
    Ib = psi(2) / Lb;
    Ic = psi(3) / Lc;
    Id = psi(4) / Ld;
    
    if Ia < 0
        Ia = 0;
    end
    if Ib < 0
        Ib = 0;
    end
    if Ic < 0
        Ic = 0;
    end
    if Id < 0
        Id = 0;
    end
    
    if Ia > 20
        Ia = 20;
    end
    if Ib > 20
        Ib = 20;
    end
    if Ic > 20
        Ic = 20;
    end
    if Id > 20
        Id = 20;
    end
    
    La = interp2(Ix,Iy,Iv,Ia,alfaMod(1),'spline'); 
    Lb = interp2(Ix,Iy,Iv,Ib,alfaMod(2),'spline');
    Lc = interp2(Ix,Iy,Iv,Ic,alfaMod(3),'spline');
    Ld = interp2(Ix,Iy,Iv,Id,alfaMod(4),'spline');
    
    Ma = interp2(Mx,My,Mv,Ia,alfaMod(1),'cubic');
    Mb = interp2(Mx,My,Mv,Ib,alfaMod(2),'cubic');
    Mc = interp2(Mx,My,Mv,Ic,alfaMod(3),'cubic');
    Md = interp2(Mx,My,Mv,Id,alfaMod(4),'cubic');
    
    Mall = Ma + Mb + Mc + Md;
    
    Mvent = 2.4691358e-07 * w^2;
    torque(i) = Mdim;
    
%     if Mall <= 0.07 && w < 10 % react
%         Mdim = 0;
%     else
%         Mdim = Mall - 0.07 - 0.05*0.05;
%     end
    
    Mdim = Mall - Mvent - 0.05*0.05;
    
%     Mdim = Mall - 0.05*0.05;
   
    w = w + Mdim / Jd * Ts;
    
    if i == 296
        i = 1 + i - 1;
    end
    
    u_ss = kss * w;
    
    speed(i) = w;
    current(i) = Id;
end
stDown = u_zs/kss - 0.02 * u_zs/kss;
(max(speed)-u_zs/kss)/(u_zs/kss)*100