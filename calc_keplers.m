function [r2,v2,tn] = calc_keplers(elements,r0,v0,xi,u,t,tol)
    a = elements(1);
    z = (xi^2)/a;
    C = calc_C(z);
    S = calc_S(z);
    tn = ((xi^3)*S+(dot(r0,v0)/sqrt(u))*(xi^2)*C+norm(r0)*xi*(1-z*S))/sqrt(u);
    if(abs(t-tn) < tol)
        tn = t;
        f = 1-(xi^2)*C/norm(r0); g = tn-(xi^3)*S/sqrt(u);
        r2 = f+g;
        fdot = sqrt(u)*xi*(z*S-1)/(norm(r0)*r2); gdot = 1-(xi^2)*C/r2;
        v2 = fdot+gdot;
        return;
    else
        r = (xi^2)*C+(dot(r0,v0)/sqrt(u))*xi*(1-z*S)+norm(r0)*(1-z*C);
        xdot = sqrt(u)/r;
        xnext = xi + (t-tn)*xdot;
        [r2,v2,tn] = calc_keplers(elements,r0,v0,xnext,u,t,tol);
        return;
    end
end