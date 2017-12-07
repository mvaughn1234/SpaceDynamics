Y = 30*pi/180
P = 10*pi/180
R = 50*pi/180

c11 = cos(P)*cos(Y);
c12 = cos(P)*sin(Y);
c13 = -sin(P);
c21 = sin(R)*sin(P)*cos(Y)-cos(R)*sin(Y);
c22 = sin(R)*sin(P)*sin(Y)+cos(R)*cos(Y);
c23 = sin(R)*cos(P);
c31 = cos(R)*sin(P)*cos(Y)+sin(R)*sin(Y);
c32 = cos(R)*sin(P)*sin(Y)-sin(R)*cos(Y);
c33 = cos(R)*cos(P);

principal_rot = acos(0.5*(c11+c22+c33-1))
principal_axis = (1/(2*sin(principal_rot)))*[c23-c32;c31-c13;c12-c21]

q1 = principal_axis(1)*sin(principal_rot/2)
q2 = principal_axis(2)*sin(principal_rot/2)
q3 = principal_axis(3)*sin(principal_rot/2)
q4 = cos(principal_rot/2)