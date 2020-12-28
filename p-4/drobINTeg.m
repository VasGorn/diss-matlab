clear

dt=1e-5;
n_point=50000;

time = (0:n_point)*dt;

mu=0.5;
c(1)=2.506628274631;
c(2)=1.000000000019;
c(3)=76.1800917294715;
c(4)=-86.5053203294168;
c(5)=24.0140982408309;
c(6)=-1.23173957245015;
c(7)=1.20865097386618E-03;
c(8)=-5.395239384953E-06;
%Riemann
x=mu;
Gamma=exp((x+0.5)*log(x+5.5)-(x+5.5)+log(c(1)*(c(2)+ c(3)/(x+1)+c(4)/(x+2)+c(5)/(x+3)+c(6)/(x+4)+c(7)/(x+5)+c(8)/(x+6))/x));
% for j=1:n_point
% k_riem(j)=(dt^x)/Gamma*(j^x-(j-1)^x)/x;
% end
for j=0:n_point
k_riem1(j+1)=(dt^x)/Gamma*((j+1)^x-(j)^x)/x;
end

for i=0:n_point
     u_in(i+1)=1;
%     if i>5000
%         u_in(i+1)=2;
%     end
    
    s=0;
    c_i=length(u_in)-1;
    for j=0:c_i
        s=s+u_in(i-j+1)*k_riem1(j+1);
    end
    y1(i+1)=s;
end

% for j=1:n_point
% k_riem2(j)=(dt^x)/Gamma*((j)^x-(j-1)^x)/x;
% end

% for i=1:n_point
%      u_in(i)=1;
% %     if i>5000
% %         u_in(i+1)=2;
% %     end
%     
%     s=0;
%     for j=1:i
%         s=s+u_in(i-j+1)*k_riem1(j);
%     end
%     y2(i)=s;
% end
%%
clear

dt=0.001;
n_point=10000;
mu=0.45;
c(1)=2.506628274631;
c(2)=1.000000000019;
c(3)=76.1800917294715;
c(4)=-86.5053203294168;
c(5)=24.0140982408309;
c(6)=-1.23173957245015;
c(7)=1.20865097386618E-03;
c(8)=-5.395239384953E-06;
%Riemann
x=mu;
Gamma=exp((x+0.5)*log(x+5.5)-(x+5.5)+log(c(1)*(c(2)+ c(3)/(x+1)+c(4)/(x+2)+c(5)/(x+3)+c(6)/(x+4)+c(7)/(x+5)+c(8)/(x+6))/x));
% for j=1:n_point
% k_riem(j)=(dt^x)/Gamma*(j^x-(j-1)^x)/x;
% end
for j=0:n_point
k_riem1(j+1)=(dt^x)/Gamma*((j+1)^x-(j)^x)/x;
end



for i=0:n_point
     u_in(i+1)=1;
    if i>5000
        u_in(i+1)=2;
    end
    
    s=0;
    c_i=length(u_in)-1;
    for j=0:c_i
        s=s+u_in(i-j+1)*k_riem1(j+1);
    end
    y1(i+1)=s;
end