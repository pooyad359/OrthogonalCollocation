function rts=legendreroots(n)
% returns roots of Legendre polynomial of degree N
p=legendrepoly(n);
rts=roots(p);