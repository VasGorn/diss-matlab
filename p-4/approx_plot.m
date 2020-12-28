%for reference w=368
numMin = 2779;
numMax = 3127;
subArray = numArray(numMin:numMax,:);
dt = 0.04;
%figure
%inertial third order
% plot((24+(0:348))*dt, (subArray(:,2)+subArray(:,3))/2)
% hold on
% grid on
% xlim([0 10])
% title('inertial third order')
% xlabel('time, s')
% plot((0:1000)*0.01,model1)

%  figure
% % %inertial second order
% plot((18+(0:348))*dt, (subArray(:,2)+subArray(:,3))/2)
% hold on
% grid on
% xlim([0 10])
% title('inertial second order')
% xlabel('time, s')
% plot((0:1000)*0.01,model2)

% figure
% %inertial second order with integrator
% plot((12+(0:348))*dt, (subArray(:,2)+subArray(:,3))/2)
% hold on
% grid on
% xlim([0 10])
% title('inertial second order with integrator')
% xlabel('time, s')
% %plot((0:1000)*0.01,model3)



%%
%for reference w=2486
dt = 0.04;
numMin = 1130;
numMax = 1255;
subArray = numArray(numMin:numMax,:);

% figure
% %inertial second order
% plot((27+(0:125))*dt, (subArray(:,2)+subArray(:,3))/2)
% hold on
% grid on
% xlim([0 5])
% title('inertial third order')
% xlabel('time, s')
% plot((0:500)*0.01,model4)

figure
%inertial second order
plot((24+(0:125))*dt, (subArray(:,2)+subArray(:,3))/2)
hold on
grid on
xlim([0 5])
title('inertial second order')
xlabel('time, s')
plot((0:500)*0.01,model5)