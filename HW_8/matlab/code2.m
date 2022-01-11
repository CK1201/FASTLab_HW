clc,clear
close all
%% 
dt=0.2;
N=20;
H=20;
R=20;
V=-0.1;
W=0.1;
path=getPath(dt,H,R,V,W);
% plot3(path(:,1),path(:,2),path(:,3))

p0=[0,0,H];
v0=[0,0,0];
a0=[0,0,0];
j=[0,0,0];

v_bound=[-6,6;-6,6;-1,6];
a_bound=[-3,3;-3,3;-1,3];
j_bound=[-3,3;-3,3;-2,2];

% p_noise=normrnd(0,0.01,size(path,1),3);
% v_noise=normrnd(0,0.01,size(path,1),3);
% a_noise=normrnd(0,0.01,size(path,1),3);
p_noise=zeros(size(path,1),3);
v_noise=zeros(size(path,1),3);
a_noise=zeros(size(path,1),3);

log=zeros(size(path,1)+1,13);
log(1,:)=[0,p0,v0,a0,j];
w1=10;
w2=1;
w3=1;
w4=1;
w5=1e4;

% Pr=path();

for i=1:size(path,1)
    t=i*dt;
    for k=1:3
%         tic
        [Tp,Tv,Ta,Bp,Bv,Ba]=getPredictionMatrix(N,dt,p0(k),v0(k),a0(k));
        H=w1*(Tp'*Tp)+w2*(Tv'*Tv)+w3*(Ta'*Ta)+w4*eye(N);
        H=blkdiag(H,w5*eye(2*N));
        if i<size(path,1)-N+1
            Bp=Bp-path(i:i+N-1,k);
        else
            Bp=Bp-[path(i:end,k);ones(N-size(path,1)+i-1,1)*path(end,k)];
        end
        F=w1*Bp'*Tp+w2*Bv'*Tv+w3*Ba'*Ta;
        F=[F,zeros(1,2*N)];
        A=[Tv,zeros(N),-eye(N);-Tv,-eye(N),zeros(N);Ta,zeros(N),zeros(N);-Ta,zeros(N),zeros(N);zeros(size(Ta)),-eye(N),zeros(N);zeros(size(Ta)),zeros(N),-eye(N);eye(N),zeros(N),zeros(N);-eye(N),zeros(N),zeros(N)];
        b=[v_bound(k,2)*ones(N,1)-Bv;-v_bound(k,1)*ones(N,1)+Bv;a_bound(k,2)*ones(N,1)-Ba;-a_bound(k,1)*ones(N,1)+Ba;zeros(N,1);zeros(N,1);j_bound(k,2)*ones(N,1);-j_bound(k,1)*ones(N,1)];
        J=quadprog(H,F,A,b);
        j(k)=J(1);
        % 模拟系统
        p0(k)=p0(k)+v0(k)*dt+0.5*a0(k)*dt^2+1/6*j(k)*dt^3+p_noise(i,k);
        v0(k)=v0(k)+a0(k)*dt+0.5*j(k)*dt^2+v_noise(i,k);
        a0(k)=a0(k)+j(k)*dt+a_noise(i,k);
        
%         toc
    end
    log(i+1,:)=[t,p0,v0,a0,j];
end
figure()
subplot(3,3,1)
plot3(path(:,1),path(:,2),path(:,3))
hold on
grid on
plot3(log(:,2),log(:,3),log(:,4),'o')
legend('参考轨迹','实际轨迹')

subplot(3,3,2)
plot(log(:,1),log(:,5))
legend('x轴速度')

subplot(3,3,3)
plot(log(:,1),log(:,7))
legend('z轴速度')

subplot(3,3,5)
plot(log(:,1),log(:,8))
legend('x轴加速度')

subplot(3,3,6)
plot(log(:,1),log(:,10))
legend('z轴加速度')

subplot(3,3,8)
plot(log(:,1),log(:,11))
legend('x轴加加速度')

subplot(3,3,9)
plot(log(:,1),log(:,13))
legend('z轴加加速度')