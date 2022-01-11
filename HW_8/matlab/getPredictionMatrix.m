function [Tp,Tv,Ta,Bp,Bv,Ba]=getPredictionMatrix(N,dt,p0,v0,a0)

Tp=zeros(N);
Tv=zeros(N);
Ta=zeros(N);

for i=1:N
    Ta(i,1:i)=ones(1,i)*dt;
end

for i=1:N
    for j=1:i
        Tv(i,j)=(i-j+0.5)*dt^2;
    end
end

for i=1:N
    for j=1:i
        Tp(i,j)=((i-j+1)*(i-j)/2+1/6)*dt^3;
    end
end

Bp=ones(N,1)*p0;
Bv=ones(N,1)*v0;
Ba=ones(N,1)*a0;

for i=1:N
    Bp(i)=Bp(i)+i*dt*v0+i^2/2*a0*dt^2;
    Bv(i)=Bv(i)+i*dt*a0;
end