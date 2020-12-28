clear
Rd = 0.68;
Uv = 4;

angleON = 7.5;
angleOFF = 22.5;

Jd = 0.0000073*2;

Ts = 1e-5;
n_point = 7000;

load('D:\cloud\DISSERT\p - 3 - data\matlab\Current=f(L,Q).mat');
load('D:\cloud\DISSERT\p - 3 - data\matlab\Torque=f(i,Q).mat');

Ix(size(Iv,1),size(Iv,2))=0;
Iy(size(Iv,1),size(Iv,2))=0;
for i = 1:size(Iv,1)
    for j = 1:size(Iv,2)
        Ix(i,j) = (j - 1) * 0.5;
        Iy(i,j) = i * 0.5;
    end
end
% i = interp2(Ix,Iy,Iv,L,Q)

Mx(size(Mv,1),size(Mv,2))=0;
My(size(Mv,1),size(Mv,2))=0;
for i = 1:size(Mv,1)
    for j = 1:size(Mv,2)
        Mx(i,j) = (j - 1) * 0.5;
        My(i,j) = (i - 1) * 0.5;
    end
end
% M = interp2(Mx,My,Mv,i,Q)

w = 0;
alfa = [0 -15 -30 -45];
alfaMod = mod(alfa, 60);
Ia = 0; Ib = 0; Ic = 0; Id = 0;
La = interp2(Ix,Iy,Iv,Ia,alfaMod(1),'spline'); 
Lb = interp2(Ix,Iy,Iv,Ib,alfaMod(2),'spline');
Lc = interp2(Ix,Iy,Iv,Ic,alfaMod(3),'spline');
Ld = interp2(Ix,Iy,Iv,Id,alfaMod(4),'spline');
psi(4) = 0.0;

Mdim = 0;

current(n_point)=0;
speed(n_point)=0;
torque(n_point)=0;

uv(n_point) = 0;
ufrac(n_point) = 0;
dE = 0.04;

time = (0:n_point-1)*Ts;
for i = 1 : n_point
    alfa = alfa + w / pi * 180 * Ts;
   
    alfaMod = mod(alfa, 60);
    sig = (alfaMod >= angleON) & (alfaMod <= angleOFF);
    u_zc = 18;
    
    sigCur = u_zc * sig;
    % current summator
%     du_c = u_zc - [Ia Ib Ic Id];
%     sig = du_c 
%     if du_c < 0.02
%         sig = (alfaMod >= angleON) & (alfaMod <= angleOFF);
%     elseif du_c > -0.02
%         sig = [0 0 0 0];
%     end

    % voltage on phase
    if (sig(1) < 1 && Ia > 0)
        Ua = -Uv;
    elseif (sig(1) > 0 && (u_zc - Ia)> dE)
        Ua = Uv;
    else
        Ua = 0;
    end
    
    if (sig(2) < 1 && Ib > 0)
        Ub = -Uv;
    elseif (sig(2) > 0 && (u_zc - Ib)> dE)
        Ub = Uv;
    else
        Ub = 0;
    end
    
    if (sig(3) < 1 && Ic > 0)
        Uc = -Uv;
    elseif (sig(3) > 0 && (u_zc - Ic)> dE)
        Uc = Uv;
    else
        Uc = 0;
    end
    
    if (sig(4) < 1 && Id > 0)
        Ud = -Uv;
    elseif (sig(4) > 0 && (u_zc - Id)> dE)
        Ud = Uv;
    else
        Ud = 0;
    end
    
    Uall = [Ua Ub Uc Ud] - [Ia Ib Ic Id] * Rd;
    psi = psi + 1 * Uall * Ts;
    
    Ia = psi(1) / La;
    Ib = psi(2) / Lb;
    Ic = psi(3) / Lc;
    Id = psi(4) / Ld;
    
    if Ia < 0
        Ia = 0;
    end
    if Ib < 0
        Ib = 0;
    end
    if Ic < 0
        Ic = 0;
    end
    if Id < 0
        Id = 0;
    end
    
    La = interp2(Ix,Iy,Iv,Ia,alfaMod(1),'spline'); 
    Lb = interp2(Ix,Iy,Iv,Ib,alfaMod(2),'spline');
    Lc = interp2(Ix,Iy,Iv,Ic,alfaMod(3),'spline');
    Ld = interp2(Ix,Iy,Iv,Id,alfaMod(4),'spline');
    
    Ma = interp2(Mx,My,Mv,Ia,alfaMod(1),'cubic');
    Mb = interp2(Mx,My,Mv,Ib,alfaMod(2),'cubic');
    Mc = interp2(Mx,My,Mv,Ic,alfaMod(3),'cubic');
    Md = interp2(Mx,My,Mv,Id,alfaMod(4),'cubic');
    
    Mall = Ma + Mb + Mc + Md;
    
    Mvent = 2.4691358e-07 * w^2;
    torque(i) = Mdim;
    
%     if Mall <= 0.07 && w < 10 % react
%         Mdim = 0;
%     else
%         Mdim = Mall - 0.07 - 0.05*0.05;
%     end
    
    Mdim = Mall - Mvent - 0.05*0.05;
    
%     Mdim = Mall - 0.05*0.05;
   
    w = w + Mdim / Jd * Ts;
    
    if i == 296
        i = 1 + i - 1;
    end
    
   
    speed(i) = w;
    current(i) = Id;
end
%(max(speed)-u_zs/kss)/(u_zs/kss)*100