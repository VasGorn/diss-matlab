load('C:\Users\kings\Documents\MATLAB\speed-loop-data\power-12v\speed2400.mat');
load('C:\Users\kings\Documents\MATLAB\speed-loop-data\power-12v\speed4000.mat');
load('C:\Users\kings\Documents\MATLAB\speed-loop-data\power-12v\speed5600.mat');
load('C:\Users\kings\Documents\MATLAB\speed-loop-data\power-12v\speed7200.mat');
load('C:\Users\kings\Documents\MATLAB\speed-loop-data\power-12v\speed8400.mat');

%%
%plot all transients
plot((1:1000-1)*0.005,speed2400(:,2)*(2*pi*200/10000));
hold on;
plot((1:1000-1)*0.005,speed4000(:,2)*(2*pi*200/10000));
plot((1:1000-1)*0.005,speed5600(:,2)*(2*pi*200/10000));
plot((1:1000-1)*0.005,speed7200(:,2)*(2*pi*200/10000));
plot((1:1000-1)*0.005,speed8400(:,2)*(2*pi*200/10000));
%arr = speed4000(1:999,2)*(2*pi*200/10000);
%csvwrite('func1.csv',arr);

