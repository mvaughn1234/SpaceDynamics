%input = [a,e,i,OMEGA,omega,nu]
function [r,v] = oe2cart(a,e,i,Omega,omega,nu)
%To be used only for planets around the sun
    input = [a,e,i,Omega,omega,nu];
    [r,v] = elementsToECI(input,132712440018);
end

function [r,v] = elementsToECI(input,u)
    i = input(3)*pi/180; OMEGA = input(4)*pi/180; omega = input(5)*pi/180;
    OM = OMEGA; om = omega;
    
    pqw = elementsToPQW(input,u);
    Rpqw = pqw(1:3);
    Vpqw = pqw(4:6);

    R11 = cos(OM)*cos(om)-sin(OM)*sin(om)*cos(i);
    R12 = -cos(OM)*sin(om)-sin(OM)*cos(om)*cos(i);
    R13 = sin(OM)*sin(i);
    R21 = sin(OM)*cos(om)+cos(OM)*sin(om)*cos(i);
    R22 = -sin(OM)*sin(om)+cos(OM)*cos(om)*cos(i);
    R23 = -cos(OM)*sin(i);
    R31 = sin(om)*sin(i);
    R32 = cos(om)*sin(i);
    R33 = cos(i);
    
    R = [[R11 R12 R13];[R21 R22 R23];[R31 R32 R33]];
    
    Rijk = R*Rpqw;
    Vijk = R*Vpqw;

    r = Rijk;
    v = Vijk;
end