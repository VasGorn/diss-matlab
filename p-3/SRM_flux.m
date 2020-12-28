clear
Rd = 0.68;
Uv = 24.2;

angleON = 7.5;
angleOFF = 22.5;

Jd = 0.0000073;

Ts = 1e-5;

n_point = 5000;

load('D:\cloud\DISSERT\p - 3 - data\matlab\Curr50.mat');
load('D:\cloud\DISSERT\p - 3 - data\matlab\Torq50.mat');

Ix(size(Iv,1),size(Iv,2))=0;
Iy(size(Iv,1),size(Iv,2))=0;
for i = 1:size(Iv,1)
    for j = 1:size(Iv,2)
        Ix(i,j) = (j - 1) * 0.0005;
        Iy(i,j) = (i - 1) * 2;
    end
end
% i = interp2(Ix,Iy,Iv,psi,Q)

Mx(size(Mv,1),size(Mv,2))=0;
My(size(Mv,1),size(Mv,2))=0;
for i = 1:size(Mv,1)
    for j = 1:size(Mv,2)
        Mx(i,j) = (j - 1) * 0.5;
        My(i,j) = (i - 1);
    end
end
% M = interp2(Mx,My,Mv,Q,i)

w = 0;
alfa = [0 -15 -30 -45];
alfaMod = mod(alfa, 60);
sig = (alfaMod >= angleON) & (alfaMod <= angleOFF);
Ia = 0; Ib = 0; Ic = 0; Id = 0;
psi(4) = 0.0;

current(n_point)=0;
speed(n_point)=0;
torque(n_point)=0;
time = (0:n_point-1)*Ts;
for i = 1 : n_point
    alfa = alfa + w / pi * 180 * Ts;
   
    alfaMod = mod(alfa, 60);
    
    sig = (alfaMod >= angleON) & (alfaMod <= angleOFF);
    
    % voltage on phase
    if (sig(1) < 1 && Ia > 0)
        Ua = -Uv;
    elseif (sig(1) > 0)
        Ua = Uv;
    else
        Ua = 0;
    end
    
    if (sig(2) < 1 && Ib > 0)
        Ub = -Uv;
    elseif (sig(2) > 0)
        Ub = Uv;
    else
        Ub = 0;
    end
    
    if (sig(3) < 1 && Ic > 0)
        Uc = -Uv;
    elseif (sig(3) > 0)
        Uc = Uv;
    else
        Uc = 0;
    end
    
    if (sig(4) < 1 && Id > 0)
        Ud = -Uv;
    elseif (sig(4) > 0)
        Ud = Uv;
    else
        Ud = 0;
    end
    
    Uall = [Ua Ub Uc Ud] - [Ia Ib Ic Id] * Rd;
    
    psi = psi + 1 * Uall * Ts;
    if psi(1) < 0
        psi(1) = 0;
    end
    if psi(2) < 0
        psi(2) = 0;
    end
    if psi(3) < 0
        psi(3) = 0;
    end
    if psi(4) < 0
        psi(4) = 0;
    end
    
    Iall = interp2(Ix,Iy,Iv,psi,alfaMod,'spline');
    if Iall(1) < 0
        Iall(1) = 0;
    end
    if Iall(2) < 0
        Iall(2) = 0;
    end
    if Iall(3) < 0
        Iall(3) = 0;
    end
    if Iall(4) < 0
        Iall(4) = 0;
    end
    
    Mall = interp2(Mx,My,Mv, alfaMod, Iall,'cubic');
    
    Mmotor = sum(Mall);
    
    Mdim = Mmotor - 0;
    
    w = w + Mdim / Jd * Ts;
    
    Ia = Iall(1); Ib = Iall(2); Ic = Iall(3); Id = Iall(4);
   
    speed(i) = w;
    current(i) = Id;
    
    if i == 138
        i = 1 + i - 1;
    end
end


