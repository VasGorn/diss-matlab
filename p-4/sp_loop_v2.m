%-------------object control--------------
%--- a1*D^(1+mu)*I + a0*D^mu*I + I = K ---
clear
a1 = 0.042015; %
a0 = 0.768569;
b0 = 1;
K = 13.555209; %
mu = 0.658944;

n_point = 1000;
dt = 0.005;
time = (0 : n_point - 1) * dt;

znam = a1*dt + a0*dt*dt;
hnK = K*dt*dt/znam;
hnA1 = (a1*dt)/znam;
hnB0 = (b0*dt*dt)/znam;

Tmu = 0.2;
kp = 24.2/10;
kds = 10/300;
%regulators
K1 = a1 / (2*Tmu*kp*kds*K);
K2 = a0 / (2*Tmu*kp*kds*K);
K3 = 1 / (2*Tmu*kp*kds*K);

%1
k_riem1(n_point) = 0;
x1=abs(mu-1);
for j=1:n_point
k_riem1(j)=(dt^x1)/gamma(x1)*(j^x1-(j-1)^x1)/x1;
end

maxVal = 12;
minVal = 0;

x(n_point) = 0;
dy(n_point) = 0;
f1(n_point) = 0;
y(n_point) = 0;
fmu1(n_point) = 0;

du_s(n_point) = 0;
u_zs(n_point) = 0;
u_1(n_point) = 0;
u_11(n_point) = 0;
u_2(n_point) = 0;
u_3(n_point) = 0;
u_sumR(n_point) = 0;
u_ds(n_point) = 0;
ua(n_point) = 0;
for i= 1:n_point
    %reference current
    u_zs(i) = 5.9524;
    %summator
    if i==1
        du_s(i) = u_zs(i)-0;
    else
        du_s(i) = u_zs(i)-u_ds(i-1);
    end
    
    %1
    s=0;
    for j=1:i
        s = s + du_s(i-j+1)*k_riem1(j);
    end
    u_11(i)= s;
    
    if i==1
         u_1(i) = K1*(u_11(i) - 0)/dt;
    else
         u_1(i) = K1*(u_11(i) - u_11(i-1))/dt;
    end
    
    if u_1(i)>maxVal %#
        u_1(i) = maxVal;
    end
    
    %2
    u_2(i)= K2*s;
%     
    if u_2(i)>maxVal
        u_2(i) = maxVal;
    end
    
    %3
    if i==1
        u_3(i)= 0 + du_s(i)*K3*dt;
    else
        u_3(i)= u_3(i-1) + du_s(i)*K3*dt;
    end
    
    if u_3(i)>maxVal %#
        u_3(i) = maxVal;
    end
    
    %reg sum
    u_sumR(i)=u_1(i)+u_2(i)+u_3(i);
    
    if u_sumR(i)>maxVal %#
        u_sumR(i) = maxVal;
    elseif u_sumR(i)<minVal
        u_sumR(i) = minVal;
    end
    
    %power converter
    if i==1
        ua(i) = 0.5*((0*Tmu + kp*u_sumR(i)*dt)/ (Tmu + dt) + (0*(Tmu - dt) + kp*u_sumR(i)*dt)/Tmu);
    else
        ua(i) = 0.5*((ua(i-1)*Tmu + kp*u_sumR(i)*dt)/ (Tmu + dt) + (ua(i-1)*(Tmu - dt) + kp*u_sumR(i)*dt)/Tmu);
    end
    %u_in(i) = kp*u_sumR(i);
    if ua(i) < 10
        x(i) = 0;
    else
        x(i) = ua(i);
    end
    
    dy(i) = hnK*x(i);
    if i > 1
        dy(i) = dy(i) + hnA1*dy(i-1) - hnB0*y(i-1);
    end
    
    f1(i) = dy(i);
    
    fmu1(i) = f1(i)*dt^mu/gamma(1+mu);
    Imu1 = 0;
    for j=1:i
        Imu1 = Imu1 + fmu1(i-j+1)*(j^mu -(j-1)^mu);
    end
    
%     if abs(Imu1) > 10
%         Imu1 = 10*sign(Imu1);
%     end
    
    y(i) = Imu1;
    
    %feedback current
    u_ds(i)= kds *Imu1;
end