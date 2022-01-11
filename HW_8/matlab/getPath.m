function path=getPath(dt,H,R,V,W)
% H=20;
% R=10;
% V=-0.2;

T=H/abs(V);
t=0:dt:T;
r=R/T*t;
x=r.*cos(W*t);
y=r.*sin(W*t);
z=H+t*V;
path=[x',y',z'];
path=[path;ones(60,3).*[x(end),y(end),z(end)]];