%-------------object control--------------
%--- a1*D^(1+mu)*I + a0*D^mu*I + I = K ---
clear
% a1 = 0.0090267; % 0.0040267
% a0 = 0.1024334;
% b0 = 1;
% K = 1.1102*4.1/24.4; % 1.1452
% mu = 0.2499;

b0 = 1;

a0 = 0.127090;
a1 = 0.00619252;
K = 0.1927788;
mu = 0.353273;

dt = 0.0002; %0.00035;%0.000008

znam = a1*dt + a0*dt*dt;
hnK = K*dt*dt/znam;
hnA1 = (a1*dt)/znam;
hnB0 = (b0*dt*dt)/znam;

n_point = 1000;

Tmu = 0.001;
kp = 24.4/10;
kdt = 10/4.1;
%regulators

% K1 = 1;% slow
% K2 = 1.59882;
% K3 = 51.9034;

K1 = 5;% fast
K2 = 1.09229;
K3 = 56.51698;

%1
%0.77553 slow
%0.794129 fast
x1=abs(0.794129);
for j=1:n_point
k_riem1(j)=(dt^x1)/gamma(x1)*(j^x1-(j-1)^x1)/x1;
end

perf(n_point) = 0;
stableTime = 0.03;
slope = 8 / stableTime;
for l = 1:n_point
    if (l - 1) * dt <= stableTime
        perf(l) = (l - 1) * dt * slope;
    else
        perf(l) = 8;
    end
end

error = 0;

x(n_point) = 0;
dy(n_point) = 0;
f1(n_point) = 0;
y(n_point) = 0;
fmu1(n_point) = 0;

du_t(n_point) = 0;
u_zt(n_point) = 0;
u_1(n_point) = 0;
u_11(n_point) = 0;
u_2(n_point) = 0;
u_3(n_point) = 0;
u_sumR(n_point) = 0;
u_dt(n_point) = 0;
u_in(n_point) = 0;

max_t = 0.0;
flag = 0;

time_in = 0.0;
for i= 1:n_point
    %reference current
    u_zt(i)= 8;
    %summator
    if i==1
        du_t(i) = u_zt(i)-0;
    else
        du_t(i) = u_zt(i)-u_dt(i-1);
    end
    
%     if i==1
%         uf(i)= (0*Tf + 1*du_t(i)*dt)/(Tf+dt);
%     else
%         uf(i)= (uf(i-1)*Tf + 1*du_t(i)*dt)/(Tf+dt);
%     end
    
    %1
    u_1(i) = K1*du_t(i);
    
%     if u_1(i)>10 %#
%         u_1(i) = 10;
%     end
    
    %2
    s=0;
    for j=1:i
        s = s + du_t(i-j+1)*k_riem1(j);
    end
    u_2(i)= K2*s;
     
%     if u_2(i)>10
%         u_2(i) = 10;
%     end
    
    %3
    if i==1
        u_3(i)= 0 + du_t(i)*K3*dt;
    else
        u_3(i)= u_3(i-1) + du_t(i)*K3*dt;
    end
    
%     if u_3(i)>10 %#
%         u_3(i) = 10;
%     end
    
    %reg sum
    u_sumR(i)=u_1(i)+u_2(i)+u_3(i);
    
    if u_sumR(i)>10 %#
        u_sumR(i) = 10;
    elseif u_sumR(i)<-10
        u_sumR(i) = -10;
    end
    
    %power converter
%     if i==1
%         u_in(i)= (0*Tmu + kp*u_sumR(i)*dt)/(Tmu+dt);
%     else
%         u_in(i)= (u_in(i-1)*Tmu + kp*u_sumR(i)*dt)/(Tmu+dt);
%     end
    if i==1
        u_in(i)= 0.5*((0*Tmu + kp * u_sumR(i)*dt)/ (Tmu + dt) + (0*(Tmu - dt) + kp * u_sumR(i)*dt)/Tmu);
    else
        u_in(i)= 0.5*((u_in(i-1)*Tmu + kp * u_sumR(i)*dt)/ (Tmu + dt) + (u_in(i-1)*(Tmu - dt) + kp * u_sumR(i)*dt)/Tmu);
    end
%     
%     u_in(i) = kp*u_in(i);
    
    dy(i) = hnK*u_in(i);
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
    %Imu1 = Imu1 + (randi(10)-5)*0.015;
    y(i) = Imu1;
    
    %feedback current
    u_dt(i)= kdt *y(i);
    
    if i > 1 && flag < 1;
        if u_dt(i)>u_dt(i-1)
            max_t = (i-1)*dt;
        else
            flag = 1;
        end
    end
    
    if i > 1
        if u_dt(i) < u_zt(1) - u_zt(1)*0.02;
            time_in = (i-1)*dt;
        end
    end
    
    error = error + (u_dt(i) - perf(i))^2;
end

error = error / n_point
ton = u_zt(1)*0.02;
up_ton = u_zt(1) + ton;
down_ton = u_zt(1) - ton;
maximum = max(u_dt);
up_ton = up_ton / kdt;
down_ton = down_ton / kdt;
maximum = maximum / kdt;
plot([0 (n_point-1)*dt],[up_ton up_ton],[0 (n_point-1)*dt],[down_ton down_ton],[0 max_t],[maximum maximum],[time_in time_in],[0 down_ton],(0:n_point-1)*dt,y)