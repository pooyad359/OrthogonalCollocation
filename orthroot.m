function r=orthroot(n)
% Calculating roots of shifted legendre polynomial (0,1)
% Input n: Quad pt rule

% Calculating the Pn(x)
% Legendre Polynomial
% Using recursive relationship
% P(order of polynomial, value of x)
% P(0,x)=1; P(1,x)=0;
% (i+1)*P(i+1,x)=(2*i+1)*x*P(i,x)-i*P(i-1,x)
syms x
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
r=double(quadpts);
r=r/2+.5;