Ff=[0.0001,0.29,0.51,0.7,0.82,0.915,0.99,1.05,1.1,1.15,1.2,1.24,1.28];
If=[0.0001,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4];
%%
%DP-12, D-12
clear
Pn=2500;
nn=1100;
Un=220;
In=16;
Rj=1.43;
N=990;
a=1;
p=2;
Rv=0.59;
wv=83;
Fn=0.0052;
nmax=3600;
Jd=0.05;

%Ff=[0.0001,0.29,0.51,0.7,0.82,0.915,0.99,1.05,1.1,1.15,1.2,1.24,1.28].*Fn*2*p*wv;
Ff=[0.0001,0.29,0.51,0.7,0.82,0.915,0.99,1.05,1.1,1.15,1.2,1.24,1.28].*Fn;
If=[0.0001,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4].*16;

Ff1=[0.0001,0.08,0.152,0.215,0.27,0.315,0.355,0.39,0.42,0.445,0.47,0.49,0.51,0.525,0.54,0.551,0.56,0.57,0.579,0.587,0.595,0.603,0.61,0.617,0.623,0.627,0.632,0.637,0.64,0.644]*Fn;%*wv;
If1=[0.0001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9]*In;


kk=p*N/2/pi/a;
Mn=kk*Fn*In;
omega_n=2*pi*nn/60;
La=0.25*Un/In*60/2/p/pi/nn;
Rsum=(Rj+Rv)*(1+0.004*(180-20));
Ta=La/Rsum;
Tv=Fn*2*p*wv/16;
Jsum=4*Jd;
Tm=Jsum*Rsum/(kk*Fn);

kp=Un/10;
Tmu=0.001;
Kdt=10/(In*3);
Kds=10/(2*pi*nmax/60);

%Tsum = Ta+Tv;
Tsum = 0.0045;

krtp1=Rsum*(Ta+Tv)/2/Tmu/kp/Kdt;
krti1=Rsum/(2*Tmu*kp*Kdt);

krti2=Rsum/(2*Tmu*Tm*kp*Kdt);

krs=Jsum*Kdt/(4*Tmu*Kds);

dt = 0.5e-4;

n_point = 40000;

ua(n_point) = 0;
de(n_point) = 0;
dpsi(n_point) = 0;
ia(n_point) = 0;
fi(n_point) = 0;
psi(n_point) = 0;
Ms = Mn;

for i=1:n_point
    ua(i) = 53.008;
        
%     if i==1
%         de(i) = ua(i)-0;
%     else
%         de(i) = ua(i)-e(i-1);
%     end
    
    if i==1
        dpsi(i) = ua(i)-0;
    else
        dpsi(i) = ua(i)-psi(i-1);
    end
    
    if i==1
        ia(i)= 0.5*((0*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (0*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
    else
        ia(i)= 0.5*((ia(i-1)*Tsum + 1 / Rsum * dpsi(i)*dt)/ (Tsum + dt) + (ia(i-1)*(Tsum - dt) + 1 / Rsum * dpsi(i)*dt)/Tsum);
    end
    
    fi(i) = spline(If1,Ff1,ia(i));
    
    if i==1
         psi(i) = 4*p*wv*(fi(i) - 0) / dt;
    else
         psi(i) = 4*p*wv*(fi(i) - fi(i-1)) / dt;
    end
    
    kF = fi(i) * kk;
    
    M(i) = ia(i) * kF *0 ;%; TEST
    
    if i==1
         Mv = 0.0019759 * 0;
    else
         Mv = 0.0019759 * w(i-1)^2;
    end
    
    Mdin = M(i) - Mv * 0;
    %Mdin = M(i) - Ms;
    
    if i==1
        w(i)= 0 + Mdin / Jsum * dt;
    else
        w(i)= w(i-1) + Mdin / Jsum * dt;
    end
    
    e(i) = w(i) * kF;
    
    
    if i==1
        iexp(i)= 0.5*((0*0.036 + 1 / Rsum * dpsi(i)*dt)/ (0.036 + dt) + (0*(0.036 - dt) + 1 / Rsum * dpsi(i)*dt)/0.036);
    else
        iexp(i)= 0.5*((iexp(i-1)*0.036 + 1 / Rsum * dpsi(i)*dt)/ (0.036 + dt) + (iexp(i-1)*(0.036 - dt) + 1 / Rsum * dpsi(i)*dt)/0.036);
    end
    
end

