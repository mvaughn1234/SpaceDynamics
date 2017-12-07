function set = orbital_elements(r,v,u)
h = cross(r,v);
n = cross([0;0;1],h);
e = (1/u)*((norm(v)^2-u/norm(r))*r-(dot(r,v))*v);
p = (norm(h)^2)/u;
a = p/(1-norm(e)^2);
OMEGA = acos(dot(n,[1;0;0])/norm(n))*180/pi;
if(n(2) < 0)
    OMEGA = - OMEGA;
end
omega = acos(dot(n,e)/(norm(n)*norm(e)))*180/pi;
if(e(3) < 0)
    omega = -omega;
end
i = acos(dot(h,[0;0;1])/norm(h))*180/pi;
nu = acos(dot(r,e)/(norm(r)*norm(e)))*180/pi;
if(dot(r,v) < 0)
    nu = -nu;
end
set = [a,norm(e),i,OMEGA,omega,nu];
end