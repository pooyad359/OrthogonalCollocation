function [x,y,p]=orthcol(fun,lbfun,rbfun,n)
%orthcol Solution of boundary condition ordinary differential equations. 
%   [x,y,p]=orthcol(fun,lbfun,rbfun,n)
%   
%   accepts the main ODE, bounadry conditions and number of orthogonal
%   points as inputs and returns the solution in form of x and y value for
%   the specific points and a polynomial that can be used to find values at 
%   the rest of the points. 
%
%   Inputs Format:
%
%   fun(x,y,dy,d2y): the ODE function
%   lbfun(x,y,dy,d2y): left boundary function
%   rbfun(x,y,dy,d2y): right boundary function
%   n: a scaler for the number of orthogonal points(length of x and y will
%   be n+2)
%   x: root of shifted Gauss-Legendre quadrature
%   y: value of the solution at each point of x
%   p: polynomial describing the solution. use polyval(p,x_new) to find the
%   value of solution at other points (x_new).
%
%   %EXAMPLE: Solving {y"-y=0 , y'(0)=0, y(1)=1} for 5 points
%
%   ode=@(x,y,dy,d2y) d2y-y
%   funl=@(x,y,dy,d2y) dy
%   funr=@(x,y,dy,d2y) y-1
%   [x,y,p]=orthcol(ode,funl,funr,5);



%%
x=fastorthroot(n);
%%

x=[0; x; 1];
n=n+2;
q=ones(n);
a=zeros(n);
b=zeros(n);

q(:,2)=x;
a(:,2)=1;
for i=3:n
    q(:,i)=x.^(i-1);
    a(:,i)=(i-1).*x.^(i-2);
    b(:,i)=(i-1)*(i-2).*x.^(i-3);
end
a=a/q;
b=b/q;
op=optimset('MaxFunEvals',1000*n,'MaxIter',1000*n);
f=fsolve(@odefun,ones(n,1),op);
aa=q\f;
p=flip(aa');
y=f;
    function r=odefun(y)
        r=fun(x,y,a*y,b*y);
        r(1)=lbfun(0,y(1),a(1,:)*y,b(1,:)*y);
        r(n)=rbfun(1,y(n),a(n,:)*y,b(n,:)*y);
    end
end