function [v1,v2] = Lamberts_Solver(r1,r2,dt,u,tm)
tolerance = 1*10^-6;
dnu = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
A = tm*sqrt(norm(r1)*norm(r2)*(1+cos(dnu)));
psi_n = 0.0;
c2 = 1/2;
c3 = 1/6;
psi_up = 4*pi^2;
psi_low = -8*pi;
check1 = true;
while check1
    check2 = true;
    yn = norm(r1) + norm(r2) + A*(psi_n*c3-1)/sqrt(c2);
    xn = sqrt(yn/c2);
    dtn = ((xn^3)*c3+A*sqrt(yn))/sqrt(u);
    if(dtn<=dt)
        psi_low = psi_n;
    else
        psi_up = psi_n;
    end
    psi_next = (psi_up+psi_low)/2;
    if(psi_next > tolerance)
        c2 = (1-cos(sqrt(psi_next)))/psi_next;
        c3 = (sqrt(psi_next)-sin(sqrt(psi_next)))/sqrt(psi_next^3);
    elseif(psi_next < -tolerance)
        c2 = (1-cosh(sqrt(-psi_next)))/psi_next;
        c3 = (sinh(sqrt(-psi_next))-sqrt(-psi_next))/sqrt((-psi_next)^3);
    else
        c2 = 1/2;
        c3 = 1/6;
    end
    psi_n = psi_next;
    if(abs(dtn-dt)<tolerance)
        check1 = false;
    end
    if(psi_low == psi_n == psi_next == psi_up)
        psi_low = psi_low*1.5;
    end
end
f = 1-yn/norm(r1);
gdot = 1-yn/norm(r2);
g = A*sqrt(yn/u);
v1 = (r2-f*r1)/g;
v2 = (gdot*r2-r1)/g;
end