function example
% EXAMPLE: Solving {y"-y=0 , y'(0)=0, y(1)=1} for 10 points
ode=@(x,y,dy,d2y) d2y-y;
funl=@(x,y,dy,d2y) dy;
funr=@(x,y,dy,d2y) y-1;
n=10;
[x,y,p]=orthcol(ode,funl,funr,n);
plot(x,y,'o-')
display(polyval(p,x))