%-------------object control--------------
%--- a1*D^(1+mu)*I + a0*D^mu*I + I = K ---
% clear
% -- 8400 --
% a1 = 0.013695; %
% a0 = 0.739272;
% b0 = 1;
% k = 14.286760; %
% mu = 0.687096;

% -- 7200 --
% a1 = 0.009788; %
% a0 = 0.76173;
% b0 = 1;
% k = 14.33070; %
% mu = 0.664222;

% -- 5600 --
% a1 = 0.042015; %
% a0 = 0.768569;
% b0 = 1;
% k =  13.555209; %
% mu = 0.658944;

% -- 4000 --
% a1 = 0.074241; %
% a0 = 0.827162;
% b0 = 1;
% k =  11.478666; %
% mu = 0.682991;

% -- 16 A --
% a1 = 7.5985838e-4; %
% a0 =  0.05246;
% b0 = 1;
% k =  3.49535; %
% mu = 0.539853;

a1 = 0.0010559; %
a0 =  0.0544679;
b0 = 1;
k =  3.438640; %
mu = 0.49718;

dt = 0.00005; %0.00035;%0.000008%655

znam = a1*dt + a0*dt*dt;
hnK = k*dt*dt/znam;
hnA1 = (a1*dt)/znam;
hnB0 = (b0*dt*dt)/znam;

n_point = 3001;
time = (0:n_point-1)*dt;

x(n_point) = 0;
dy(n_point) = 0;
f1(n_point) = 0;
y(n_point) = 0;
fmu1(n_point) = 0;

sum = 0;

for i= 1:n_point
    x(i)= 5.0;
    %x(i)= 20.74;
%     x(i)= 16.133;
%     x(i)= 11.52;
    
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
    
%     er = (y(i) - y2(i))^2;
%     sum = sum + er;
end
