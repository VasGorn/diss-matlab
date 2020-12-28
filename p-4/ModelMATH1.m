%------object control-----
%--- a0*D^mu*I + I = K ---
clear
% -- 8400 --
% a1 = 0;
% a0 = 0.72823;
% b0 = 1;
% k = 14.174321;
% mu = 0.692769;

% -- 7200 --
% a1 = 0;
% a0 = 0.731159;
% b0 = 1;
% k = 13.90429;
% mu = 0.712554;

% -- 5600 --
% a1 = 0;
% a0 = 0.709188;
% b0 = 1;
% k = 13.02072;
% mu = 0.73234;

% -- 4000 --
% a1 = 0;
% a0 = 0.681359;
% b0 = 1;
% k = 10.42127;
% mu = 0.81954;

% -- 2400 --
% a1 = 0;
% a0 = 0.49387;
% b0 = 1;
% k = 3.64094;
% mu = 0.95070;

% -- 16 A --
a1 = 0;
a0 = 0.003992;
b0 = 1;
k = 2.998968;
mu = 1.21632;

dt = 0.00005;

znam = a1*dt + a0*dt*dt;
hnK = k*dt*dt/znam;
hnA1 = (a1*dt)/znam;
hnB0 = (b0*dt*dt)/znam;

n_point = 2001;
time = (0:n_point-1)*dt;

x(n_point) = 0;
dy(n_point) = 0;
f1(n_point) = 0;
y(n_point) = 0;
fmu1(n_point) = 0;
for i= 1:n_point
%     x(i)= 24.2;
%     x(i)= 20.74;
%     x(i)= 16.133;
%     x(i)= 11.52;
    x(i)= 5;
    
    dy(i) = hnK*x(i);
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
    
end