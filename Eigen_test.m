close all
clear
clc

v1 = [0 0 1]';
% v1 = rand(3,1);
% v1 = v1/norm(v1);
v2 = [sqrt(3)/2 1/2 0]';
p1 = [1 0 0]';
p2 = [-1/2 sqrt(3)/2, 0]';
r(:,1) = [v1; 0];
r(:,2) = [v2; 0];
r(:,3) = [p1; 1];
r(:,4) = [p2; 1];



k  = [2 2 2 2];% abs(rand(1,4));%


[~,n] = size(r);
AA  = zeros(4);
for i=1:n
    AA = AA + k(i)*r(:,i)*r(:,i)';
end

A = AA(1:3,1:3);
b = AA(1:3,4);
d = AA(4,4);

Q = A - (1/d)*b*b';
DD = [Q zeros(3,1);zeros(1,3) d];


WA = trace(A)*eye(3) - A;

MA = [-0.5*WA -0.5*skew(b);
      skew(b) -d*eye(3)]
  
  
WQ = trace(Q)*eye(3) - Q;
MQ = [-0.5*WQ zeros(3);
     zeros(3) -d*eye(3)]
 



t= 0:0.01:5;
u = 15*randn(length(t),6);
x0 = [randn(3,1);10*randn(3,1)];


sysA = ss(2*MA,eye(6),eye(6),zeros(1)); 
sysQ= ss(2*MQ,eye(6),eye(6),zeros(1));   
yA = lsim(sysA,u,t,x0);
yQ= lsim(sysQ,u,t,x0);

figure
subplot(2,1,1)
plot(t,sqrt(sum(yA(:,1:3).^2,2)),'k--',t,sqrt(sum(yQ(:,1:3).^2,2)),'r-.','linewidth',2)
xlabel('time (s)')
ylabel('angle')
subplot(2,1,2)
plot(t,sqrt(sum(yA(:,4:6).^2,2)),'k--',t,sqrt(sum(yQ(:,4:6).^2,2)),'r-.','linewidth',2)
xlabel('time (s)')
ylabel('position')

