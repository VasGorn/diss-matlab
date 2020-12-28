clear
Rd = 0.68;
Umax = 24.2;

angleON = 7.5;
angleOFF = 22.5;

Jd = 0.0000073 * 2;

Ts = 1e-5;
n_point = 60000;
time = (0:n_point-1)*Ts;

load('D:\cloud\DISSERT\p - 3 - data\matlab\Current=f(L,Q).mat');
load('D:\cloud\DISSERT\p - 3 - data\matlab\Torque=f(i,Q).mat');

Ix(size(Iv,1),size(Iv,2))=0;
Iy(size(Iv,1),size(Iv,2))=0;
for i = 1:size(Iv,1)
    for j = 1:size(Iv,2)
        Ix(i,j) = (j - 1) * 0.5;
        Iy(i,j) = i * 0.5;
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

% a1 = 0;
% a0 = 0.081731802;
% b0 = 1;
% k = 33.472190;
mu = 0.71;

Tmus = 0.005;
kss = 10 / 500;
Kp = 24.2 / 10;
dt = Ts;

% speed coeff
mu_c = 1 + mu;
% a = 1.1;
% b=exp(-10.27+7.831*mu_c)/a/1.5;
% a = 1.6;
% b=exp(-10.27+7.831*mu_c)/a;

% K0s = 1 / (a * b * Tmus^(1+mu) * Kp * kss * k);
% K1s = b * Tmus * a0 * K0s;
% K2s = b * Tmus * K0s;
% K3s = a0 * K0s;
% K4s = 1 * K0s*0;

k_riem1(n_point) = 0;
x1=abs(mu);
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

uv1(n_point) = 0;
uv2(n_point) = 0;
uv3(n_point) = 0;
uv4(n_point) = 0;
uv5(n_point) = 0;

ufrac(n_point) = 0;
u_c = 0;
u_ss = 0;
u3 = 0;
u2 = 0;
u4 = 0;

% x = [0 4 5 6 7 8 9 10];
% y = [1.2 1.2 1 0.9 0.85 0.8 0.7 0.7];
x = [0 4 9 10];
y = [1.2 1.2 0.7 0.7];


for i = 1 : n_point
    alfa = alfa + w / pi * 180 * Ts;
    
    u_zs = 3;
    
    if i > 20000 && i <= 40000
        u_zs = 6;
    elseif i > 40000 && i <= 60000
        u_zs = 9;
    end
    
    % speed summator
    du_s = u_zs - u_ss;
    
    % change coeff
%     a = 1.1;
%     b=exp(-10.27+7.831*mu_c)/a;
    if u_zs >= 8.33 
        a = 1.18;
        b=exp(-10.27+7.831 * mu_c) / a / 1.7;
  
        a0 = 0.057960;
        k = 25.49716;
    elseif u_zs >= 6.67
        a = 1.15;
        b=exp(-10.27+7.831 * mu_c) / a / 1.6;
  
        a0 = 0.05850;
        k = 26.9419;
    elseif u_zs >= 5
        a = 1.2;
        b=exp(-10.27+7.831 * mu_c) / a * 1.2;
        
        a0 = 0.0601538;
        k = 28.83720;
    elseif u_zs >= 3.33
        a = 1.3;
        b=exp(-10.27+7.831*mu_c)/a;
        
        a0 = 0.075079;
        k = 33.281450;    
    elseif u_zs >= 1.67
        a = 1.4;
        b=exp(-10.27+7.831*mu_c)/a*1.4;
 
        a0 = 0.1071844;
        k = 41.7249;
    elseif u_zs >= 0
        a = 1.6;
        b=exp(-10.27+7.831*mu_c)/a * 1.3;
        
        a0 = 0.29512428;
        k = 80.84416;
    end
    
    K0s = 1 / (a * b * Tmus^(1+mu) * Kp * kss * k);
    K1s = b * Tmus * a0 * K0s;
    K2s = b * Tmus * K0s;
    K3s = a0 * K0s;
    K4s = 1 * K0s / 5;
    
    % 1 
    u1 = K1s * du_s;
%     if u1>10
%         u1 = 10;
%     elseif u1<-10
%         u1 = -10;
%     end
    
    % 2
    ufrac(i) = du_s;
    uv5(i) = du_s;
    s=0;
    for j=1:i
        s = s + ufrac(i-j+1)*k_riem1(j);
    end
    u2 = s * K2s;

    % 3
    u3 = u3 + K3s * du_s * dt;
    
    % 4
    u4 = u4 + K4s * s * dt;
    if u4>10
        u4 = 10;
    elseif u4<-10
        u4 = -10;
    end
    
    uv1(i) = u1;
    uv2(i) = u2;
    uv3(i) = u3;
    uv4(i) = u4;
    
    ureg = u1 + u2 + u3 + u4;
    
    u_c= 0.5*((u_c * Tmus + Kp * ureg*dt)/ (Tmus + dt) + (u_c*(Tmus - dt) + Kp * ureg*dt)/Tmus);
    
    if u_c > 0
        angleON = 7.5;
        angleOFF = 22.5;
    else
        angleON = 37.5;
        angleOFF = 52.5;
    end
    
    alfaMod = mod(alfa, 60);
    
    
    sig = (alfaMod >= angleON) & (alfaMod <= angleOFF);
    Uv = abs(u_c);
%     Uv = 6;
    
    % voltage on phase
    if (sig(1) < 1 && Ia > 0)
        Ua = -Umax;
    elseif (sig(1) > 0)
        Ua = Uv;
    else
        Ua = 0;
    end
    
    if (sig(2) < 1 && Ib > 0)
        Ub = -Umax;
    elseif (sig(2) > 0)
        Ub = Uv;
    else
        Ub = 0;
    end
    
    if (sig(3) < 1 && Ic > 0)
        Uc = -Umax;
    elseif (sig(3) > 0)
        Uc = Uv;
    else
        Uc = 0;
    end
    
    if (sig(4) < 1 && Id > 0)
        Ud = -Umax;
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
    
    if i == 111
        i = 1 + i - 1;
    end
    
    u_ss = kss * w;
    
    speed(i) = w;
    current(i) = Id;
end
stDown = u_zs/kss - 0.02 * u_zs/kss;
stStabe = u_zs/kss;
(max(speed)-u_zs/kss)/(u_zs/kss)*100