function E=findE(M,e)
%%%Find the eccentric anomaly from the mean motion and the eccentricity, M
%%%should be in radians
%%%% created by Dr. Christine Hartzell (2010)
%%%% Distributed to ENAE 404 for use in project assignment

if M>-pi&&M<0 ||M>pi
    E=M-e;
else
    E=M+e;
end

step=(M-E+e*sin(E))/(1-e*cos(E));
while abs(step)>1E-13
    E=E+step;
    step=(M-E+e*sin(E))/(1-e*cos(E));
end