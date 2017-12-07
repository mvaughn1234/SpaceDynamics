function xi = xiGuess(elements,t,u)
    a = elements(1);
    period = 2*pi*sqrt((a^3)/u);
    xi = 2*pi*sqrt(a)*(t/period);
end