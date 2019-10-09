function quadpts30
% This function only prints the roots of Gauss-Legendre quadrature (n = 0 -
% 30). ***FOR TEST ONLY.***
% Use orthroot to get the values returned.
syms x

for n=1:30
    m=n-1;
    P0=1;
    P1=x;
    for i=1:1:m
        Pn=((2.0*i+1)*x*P1-i*P0)/(i+1.0);
        P0=P1;
        P1=Pn;
    end
    if n==1
        Pn=P1;
    end
    Pn=expand(Pn);
    quadpts=solve(vpa(Pn,32));
    quadpts=sort(quadpts);
    disp('r=[')
    disp(quadpts)
    disp('];')
end
