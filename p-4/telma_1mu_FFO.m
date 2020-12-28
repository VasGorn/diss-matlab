%-------------object control--------------
%--- a1*D^(1+mu)*I + a0*D^mu*I + I = K ---
clear

a1 = 2.8951; %
a0 = 1.8987;
b0 = 1;
K =  0.11514; %
mu = 0.626101;

% for 449
% a1 = 2.790224; %
% a0 = 2.691900;
% b0 = 1;
% K =  0.078585; %
% mu = 0.608827;

n_point = 1000;
dt = 0.01;
time = (0 : n_point - 1) * dt;

znam = a1*dt + a0*dt*dt;
hnK = K*dt*dt/znam;
hnA1 = (a1*dt)/znam;
hnB0 = (b0*dt*dt)/znam;

Tf = 0.1;
kf = 1;
mu_c = mu + 1;

%regulators
a=0.45;
%a=0.25;
b=35;

% K0 = 1 / (a * Tf^(mu_c - 1));
% K1 = a1 / (1 * 1 *K);
% K2 = (a0*b*Tf + a1) / (b * Tf * 1 * 1 * K);
% K3 = a0 / (b * Tf * 1 * 1 * K);
% K4 = 1 / (b * Tf * 1 * 1 * K);
% K5 = 1 / (1 * 1 * K);

K0 = 1 / (a * Tf^(mu_c - 1));
K1 = K0 * a1 / (1 * 1 *K);
K2 = K0 * (a0*b*Tf + a1) / (b * Tf * 1 * 1 * K);
K3 = K0 * a0 / (b * Tf * 1 * 1 * K);
K4 = K0 * 1 / (b * Tf * 1 * 1 * K);
K5 = K0 * 1 / (1 * 1 * K);

%1
x1=abs(mu);
for j=1:n_point
k_riem(j)=(dt^x1)/gamma(x1)*(j^x1-(j-1)^x1)/x1;
end

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
u_4(n_point) = 0;
u_5(n_point) = 0;
u_sumR(n_point) = 0;
u_ds(n_point) = 0;
ua(n_point) = 0;
uf(n_point) = 0;

Imu1 = 0.0;
for i= 1:n_point
    %reference current
    u_zs(i) = 100;
    %summator
    if i==1
        du_s(i) = u_zs(i)-0;
    else
        du_s(i) = u_zs(i)-Imu1;
    end
    
    %1---------------------------------------------------------   
    if i==1
         u_1(i) = K1*(du_s(i) - 0) / dt;
    else
         u_1(i) = K1*(du_s(i) - du_s(i-1))/dt;
    end
    
    %2---------------------------------------------------------
    u_2(i) = K2 * du_s(i);
    
    %3---------------------------------------------------------
    if i==1
        u_3(i)= 0 + du_s(i) * K3 * dt;
    else
        u_3(i)= u_3(i-1) + du_s(i) * K3 * dt;
    end
    
%     if u_3(i)>10 %#
%         u_3(i) = 10;
%     end
    
    %4---------------------------------------------------------
    s=0;
    for j=1:i
        s = s + du_s(i-j+1)*k_riem(j);
    end
    
    if i==1
        u_4(i)= 0 + s * K4 * dt;
    else
        u_4(i)= u_4(i-1) + s * K4 * dt;
    end
    
%     if u_4(i)>10 %#
%         u_4(i) = 10;
%     end
    
    %5---------------------------------------------------------
    u_5(i) = K5 * s;
    
%     if u_5(i)>10
%         u_5(i) = 10;
%     end
    
    %reg sum
    u_sumR(i)=u_1(i) + u_2(i) + u_3(i) + u_4(i) + u_5(i);
    
    if i==1
        uf(i) = 0.5*((0*Tf + kf*u_sumR(i)*dt)/ (Tf + dt) + (0*(Tf - dt) + kf*u_sumR(i)*dt)/Tf);
    else
        uf(i) = 0.5*((uf(i-1)*Tf + kf*u_sumR(i)*dt)/ (Tf + dt) + (uf(i-1)*(Tf - dt) + kf*u_sumR(i)*dt)/Tf);
    end
    
%     if(uf(i) > 42000)
%         uf(i) = 42000;
%     end
%     if(uf(i) < 0)
%         uf(i) = 0;
%     end
    
    dy(i) = hnK*uf(i);
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
   
end

%%
% load data for Kp50_Ki50_Kd9
clear
filePath = 'C:\Users\kings\Documents\MATLAB\telma_19_12_05\';
fileName = 'Kp50_Ki50_Kd9_undep.csv';
fullFile = strcat(filePath, fileName);

delimiter = ';';
formatSpec = '%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%[^\n\r]';
fileID = fopen(fullFile,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);

numArray = [dataArray{1:end-1}];
%clearvars -except numArray

% for 368
numMin = 2779;
numMax = 3127;
subArray = numArray(numMin:numMax,:);
dt = 0.04;

%figure
%fractinal order 1 + mu
plot((15+(0:348))*dt, (subArray(:,2)+subArray(:,3))/2)
hold on
grid on
title('fractinal order 1 + mu')
xlabel('time, s')
%plot((0:1000)*0.01,model3)

%%
% load data for Kp25_Ki25_Kd1_5_undep
clear
filePath = 'C:\Users\kings\Documents\MATLAB\telma_19_12_05\';
fileName = 'Kp25_Ki25_Kd1_5_undep.csv';
fullFile = strcat(filePath, fileName);

delimiter = ';';
formatSpec = '%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%[^\n\r]';
fileID = fopen(fullFile,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);

numArray = [dataArray{1:end-1}];
%clearvars -except numArray

% for 449
numMin = 2350;
numMax = 2700;
subArray = numArray(numMin:numMax,:);
dt = 0.04;

%figure
%fractinal order 1 + mu
plot((28+(0:350))*dt, (subArray(:,2)+subArray(:,3))/2)
hold on
grid on
title('Kp=25, Ki=25, Kd=1.5, fractinal order 1 + mu, a1=2.79, a0=2.692, K=0.0786, mu=0.609')
xlabel('time, s')
