% load data
clear
filePath = 'C:\Users\kings\Documents\MATLAB\telma_19_12_25\';
fileName = 'Kd_210_Kp_248_Ki_44_Kffi_16_Kfi_75_Tf_0_2_limit_In-30_filterFront_LIN_3P.csv';
fullFile = strcat(filePath, fileName);

delimiter = ';';
% formatSpec = '%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%[^\n\r]';
%formatSpec = '%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%[^\n\r]';
% formatSpec = '%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%*q%f%[^\n\r]';
% current
formatSpec = '%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%[^\n\r]';

fileID = fopen(fullFile,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);

numArray = [dataArray{1:end-1}];

clearvars -except numArray
%% plot all data
PLOT_NAME = 'Kp = 65, Ki = 50, Kd = 15';
dt = 0.04;

YMAX = 3000;
YMIN = 0;
STEP = 100;

%STEP_TIME = 5;

n = length(numArray);
maxTime = (n - 1) * dt;
time = transpose((0:n-1)*dt);

figure('Color','w')
plot((0:n-1)*dt, numArray(:,1)*0.6, 'LineWidth',2)
hold on
grid on

% Create ylabel
ylabel('n, rpm','FontName','Times New Roman');
% Create xlabel
xlabel('time, s','FontName','Times New Roman');


plot((0:n-1)*dt, abs(numArray(:,2))*0.6, 'LineWidth',2)
plot((0:n-1)*dt, abs(numArray(:,3))*0.6, 'LineWidth',2)
plot((0:n-1)*dt, (abs(numArray(:,2)) +abs(numArray(:,3))) / 2*0.6 , 'LineWidth',2)

%ylim([YMIN YMAX])
%xlim([0 maxTime])
%set(gca,'XTick',0:STEP_TIME:maxTime);
%set(gca,'YTick',YMIN:STEP:YMAX);

legend('n_{\itref}', 'n_{\itleft}', 'n_{\itright}', 'n_{\itave}')

figure
plot(numArray(:,1))
hold on
grid on
plot(numArray(:,2))

figure('Color','w')
plot((0:n-1)*dt, numArray(:,6))
hold on
grid on

title(PLOT_NAME)
xlabel('time, s')
ylabel('Torque, kg\cdotm')
%set(gca,'XTick',35:1:47);
%xlim([35.2 47.2])
plot((0:n-1)*dt, numArray(:,7))
plot((0:n-1)*dt, ((numArray(:,6) + numArray(:,7)) / 2))
legend('T_{\itleft}', 'T_{\itright}', 'T_{\itave}')

%% plot TRANSIENTS and calculate ROOT MEAN SQUARE ERORR = MARK
numMin = 1052;
numMax = 1400;

YMIN = 2600;
YMAX = 2400;
STEP = 20;

subArray = numArray(numMin:numMax,:);
ns = length(subArray);
minTimeS = numMin * dt;
maxTimeS = numMax * dt;
time_s = transpose((numMin:numMax) * dt);

figure('Color','w')
plot((0:ns-1)*dt, 0.6*subArray(:,1))
hold on
grid on

xlabel('time, s')
ylabel('n, rpm')

% plot((0:ns-1)*dt, 0.6*abs(subArray(:,2)))
% plot((0:ns-1)*dt, 0.6*abs(subArray(:,3)))
plot((0:ns-1)*dt, 0.6*(abs(subArray(:,2)) + abs(subArray(:,3))) / 2)

%ylim([YMIN YMAX])
%xlim([minTimeS maxTimeS])
%set(gca,'XTick',0:1:maxTimeS);
%set(gca,'YTick',YMIN:STEP:YMAX);

% valSum = 0;
% for i=1 : ns
%     x1 = (subArray(i, 1) - subArray(i, 2))^2;
%     x2 = (subArray(i, 1) - subArray(i, 3))^2;
%     
%     valSum = valSum + x1 + x2;
% end

% valQ = sqrt(valSum / ns)
% 
% NAME_S = strcat(PLOT_NAME, ', MARK= ', num2str(valQ));
% title(NAME_S)

%clearvars -except numArray dt PLOT_NAME;

