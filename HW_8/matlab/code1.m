clc,clear
close all
%% 
p0=-10;
v0=-3;
a0=0;
N=20;
dt=0.2;
log=[0,p0,v0,a0];
w1=1;
w2=1;
w3=1;
w4=1;
w5=1e4;
Pr=ones(N,1)*2;
for t=0.2:0.2:20
    [Tp,Tv,Ta,Bp,Bv,Ba]=getPredictionMatrix(N,dt,p0,v0,a0);
    H=w1*(Tp'*Tp)+w2*(Tv'*Tv)+w3*(Ta'*Ta)+w4*eye(N);
    H=blkdiag(H,w5*eye(2*N));
    Bp=Bp-Pr;
    F=w1*Bp'*Tp+w2*Bv'*Tv+w3*Ba'*Ta;
    F=[F,zeros(1,2*N)];
    A=[Tv,zeros(N),-eye(N);-Tv,-eye(N),zeros(N);Ta,zeros(N),zeros(N);-Ta,zeros(N),zeros(N);zeros(size(Ta)),-eye(N),zeros(N);zeros(size(Ta)),zeros(N),-eye(N)];
    b=[2*ones(N,1)-Bv;ones(N,1)+Bv;ones(N,1)-Ba;ones(N,1)+Ba;zeros(N,1);zeros(N,1)];
    J=quadprog(H,F,A,b);
    j=J(1);
    % 模拟系统
    p0=p0+v0*dt+0.5*a0*dt^2+1/6*j*dt^3;
    v0=v0+a0*dt+0.5*j*dt^2;
%     if t==4
%         v0=v0+4;
%     end
    a0=a0+j*dt;
    log=[log;t,p0,v0,a0];
end
figure()
hold on
grid on
plot(log(:,1),log(:,2))
plot(log(:,1),log(:,3))
plot(log(:,1),log(:,4))
legend('p','v','a')