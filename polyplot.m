function polyplot(p,a,b)
t=linspace(a,b,1000);
y=polyval(p,t);
plot(t,y)