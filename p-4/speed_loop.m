clear
% control object
k = 13.021;
a0 = 0.709;
mu = 0.732;

a1 = 0;
b0 = 1;

n_point = 1000;
dt = 0.005;
time = (0:n_point-1)*dt;

znam = a1*dt + a0*dt*dt;
hnK = k*dt*dt/znam;
hnA1 = (a1*dt)/znam;
hnB0 = (b0*dt*dt)/znam;

Tmu = 0.2;
ks = 10/ 300;
kconv = 24.2 / 10;

k_riem1(n_point) = 0;
x1=abs(1 - mu);
for j=1:n_point
    k_riem1(j)=(dt^x1)/gamma(x1)*(j^x1-(j-1)^x1)/x1;
end

K1 = a0 / (2 * Tmu * k * ks * kconv);
K2 = 1 / (2 * Tmu * k * ks * kconv);

u_zs(n_point) = 0;
du_s(n_point) = 0;
u_1s(n_point) = 0;
u_2s(n_point) = 0;
u_s(n_point) = 0;
u_ds(n_point) = 0;
ua(n_point) = 0;

dy(n_point) = 0;
f1(n_point) = 0;
fmu1(n_point) = 0;
y(n_point) = 0;
x(n_point) = 0;

for i=1:n_point
    u_zs(i) = 5.9524;
    
    % speed summator
    if i==1
        du_s(i) = u_zs(i);
    else
        du_s(i) = u_zs(i)-u_ds(i-1);
    end
    
    % speed regulator
    %-1- I^(2-mu)   -----------------------------------    
    s=0;
    for j=1:i
        s = s + du_s(i-j+1)*k_riem1(j);
    end
    
    u_1s(i) = K1 * s;
    
    if u_1s(i) > 10 %#
        u_1s(i) = 10;
    elseif u_1s(i) < 0
        u_1s(i) = 0;
    end
    
    %-2- I     -----------------------------------
    if i==1
        u_2s(i)= du_s(i)*K2*dt;
    else
        u_2s(i)= u_2s(i-1) + du_s(i)*K2*dt;
    end
    
    
    if u_2s(i) > 10 %#
        u_2s(i) = 10;
    elseif u_2s(i) < 0
        u_2s(i) = 0;
    end
    
    u_s(i) = u_1s(i) + u_2s(i);
    
    if u_s(i) > 10 %#
        u_s(i) = 10;
    elseif u_s(i) < 0
        u_s(i) = 0;
    end
    
    % power converter    
    if i==1
        ua(i) = 0.5*((0*Tmu + kconv*u_s(i)*dt)/ (Tmu + dt) + (0*(Tmu - dt) + kconv*u_s(i)*dt)/Tmu);
    else
        ua(i) = 0.5*((ua(i-1)*Tmu + kconv*u_s(i)*dt)/ (Tmu + dt) + (ua(i-1)*(Tmu - dt) + kconv*u_s(i)*dt)/Tmu);
    end
    if ua(i) < 10
        x(i) = 0;
    else
        x(i) = ua(i);
    end
    
    dy(i) = hnK * x(i);
    if i > 1
        dy(i) = dy(i) + hnA1*dy(i-1) - hnB0*y(i-1);
    end
    
    f1(i) = dy(i);
    
    fmu1(i) = f1(i)*(dt^mu)/gamma(1+mu);
    Imu1 = 0;
    for j=1:i
        Imu1 = Imu1 + fmu1(i-j+1)*(j^mu -(j-1)^mu);
    end
    
%     if abs(Imu1) > 10
%         Imu1 = 10*sign(Imu1);
%     end
    
    y(i) = Imu1;
    
    u_ds(i) = Imu1 * ks;
    
end
