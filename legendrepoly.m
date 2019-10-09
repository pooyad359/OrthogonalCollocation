function p=legendrepoly(n)
% Produces Legendre polynomial of degree N

if n==0
    p=1;
elseif n==1
    p=[1,0];
elseif n==2
    p=[3/2,0,-1/2];
elseif n==3
    p=[2.5,0,-1.5,0];
elseif n==4
    p=[4.3750, 0 ,-3.7500,0, 0.3750];
elseif n>4
    p0=legendrepoly(n-2);
    p1=legendrepoly(n-1);
    p0=-(n-1)/n*p0;
    p1=conv([(2*n-1)/n 0],p1);
    p=polysum(p0,p1);
else
    error('invalid input')
end

function p=polysum(p1,p2)
n1=length(p1);
n2=length(p2);
if n1==n2
    p=p1+p2;
    return
end
if n1>n2
    pb=p1;
    ps=p2;
else
    pb=p2;
    ps=p1;
end
dn=abs(n1-n2);
p=pb+[zeros(1,dn) ps];
    