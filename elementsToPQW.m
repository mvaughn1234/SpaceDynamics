%input = [a,e,i,OMEGA,omega,nu]

function PQW = elementsToPQW(input,u)
    a = input(1); e = input(2); i = input(3)*pi/180;
    OMEGA = input(4)*pi/180; omega = input(5)*pi/180; nu = input(6)*pi/180;
    
    p = a*(1-e^2);
    
    rnorm = p/(1+e*cos(nu));
    r = [rnorm*cos(nu);rnorm*sin(nu);0];

    v = sqrt(u/p)*[-sin(nu);(e+cos(nu));0];
    
    PQW = [r;v];
end

