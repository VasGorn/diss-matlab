%------object control-----
%--- a0*D^mu*I + I = K ---
%clear
% 24 V
% uv = 24;
% a1 = 0;
% a0 = 0.036638;
% b0 = 1;
% k = 23.57993;
% mu = 0.78489;

% 20 V
% uv = 20;
% a1 = 0;
% a0 = 0.045715;
% b0 = 1;
% k = 25.63896;
% mu = 0.739205;

% 16 V
% uv = 16;
% a1 = 0;
% a0 = 0.093781;
% b0 = 1;
% k = 31.44662;
% mu = 0.62896;

% 12 V
% uv = 12;
% a1 = 0;
% a0 = 0.104321;
% b0 = 1;
% k = 36.69840;
% mu = 0.66019;

% 8 V
% uv = 8;
% a1 = 0;
% a0 = 0.11412;
% b0 = 1;
% k = 42.86211;
% mu = 0.70981;

% 4 V
uv = 4;
a1 = 0;
a0 = 0.0487007;
b0 = 1;
k = 40.81853;
mu = 0.95482;

%===========================================
% mu = 0.71;
% b0 = 1;
% a1 = 0;

% 24 V
% uv = 24;
% a0 = 0.057960;
% k = 25.49716;

% 20 V
% uv = 20;
% a0 = 0.05850;
% k = 26.9419;

% 16 V
% uv = 16;
% a0 = 0.0601538;
% k = 28.83720;

% 12 V
% uv = 12;
% a0 = 0.075079;
% k = 33.281450;

% 8 V
% uv = 8;
% a0 = 0.1071844;
% k = 41.7249;

% 4 V
% uv = 4;
% a0 = 0.29512428;
% k = 80.84416;

dt = 1e-5; %0.00035;%0.000008%655

znam = a1*dt + a0*dt*dt;
hnK = k*dt*dt/znam;
hnA1 = (a1*dt)/znam;
hnB0 = (b0*dt*dt)/znam;

n_point = 7000;
p=[-849254294.075150,1080882839.78999,-583864324.415683,173757781.706893,-30851093.9552781,3282760.98286081,-194590.447964767,4727.72312141719,61.0492854865135,0.131889338497552];
poly(n_point) = 0;

for l = 1:n_point
    x = (l-1)*dt;
    poly(l) = polyval(p,x);
    
end

error = 0;

x(n_point) = 0;
dy(n_point) = 0;
f1(n_point) = 0;
y(n_point) = 0;
fmu1(n_point) = 0;

time = (0:n_point-1)*dt;

for i= 1:n_point
    x(i)= uv;
    
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
    %if(i > 700)
        %error = error + (speed(i)/500 - Imu1/500)^2;
    %end
    
end
%plot((0:n_point-1)*dt,y,(0:n_point-1)*dt, polyval(p,(0:n_point-1)*dt))
% e = error
% h = (error / (n_point))
% g = (error / (n_point))^0.5