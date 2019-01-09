
clear all
close all
clc

global r k ba index AA DD U theta_q delta delta_m k_beta 

%% setup
v1 = [0 0 1]';
v2 = [sqrt(3)/2 1/2 0]';
p1 = 2*[0 0 1]';
p2 = p1+[-1/2, sqrt(3)/2, 0]';
p3 = p1+[sqrt(2)/2, sqrt(2)/2, 0]';
%simulation 1
r(:,1) = [v1; 0];
r(:,2) = [v2; 0];
r(:,3) = [p3; 1];
v3 = cross(v1,v2);
v3 = v3/norm(v3);
r(:,4) = [v3; 0];

% % simulation 2
% r(:,1) = [p1; 1];
% r(:,2) = [p2; 1];
% r(:,3) = [p3; 1];
% pc = (p1+p2+3)/3;
% v = cross(p1-pc,p2-pc);
% v = v/norm(v);
% r(:,4) = [v; 0];


k  = [1 1 1 1];



theta_q = 2*pi/3;
delta = 2*(1-cos(theta_q))*0.5; %0.5*4*(1-cos(theta_q))/3;  % 

k_beta = 1;
% T= [0.8*eye(3) zeros(3);zeros(3) 40*eye(3)];

bw    = [-0.02 0.02 0.1]';
bv    = [0.02 -0.1 0.01]';
ba = [bw; bv];

index = 0;
[~,n] = size(r);
AA  = zeros(4);
for i=1:n
    AA = AA + k(i)*r(:,i)*r(:,i)';
end

A = AA(1:3,1:3);
B = AA(1:3,4);
d = AA(4,4);

QQ = A - (1/d)*B*B';
DD = [QQ zeros(3,1);zeros(1,3) d];

[U,~] = eig(QQ);






%%  initial conditions
% Initialization for real trajectory
% u      = [1 0 -1]';
% u = randn(3,1);
u      =   U(:,1);% u/norm(u); % 
theta  = pi/10;
Q_0  = [cos(theta/2) sin(theta/2)*u']';
p_0  = [0 0 15]'; % (eye(3)-quat2rotm(Q_0'))*(1/d)*B; % 

% Initialization for observer
Qhat_0 = [1 0 0 0]';
phat_0  =  [0 0 0]'; 
bahat_0 = [zeros(3,1);0;0;0];
x0 = [Q_0;p_0;Qhat_0;phat_0;bahat_0;0];


delta_m = 0.00;

% simulation horizon
TSPAN=[0 100];
JSPAN = [0 5];
% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;
options = odeset('RelTol',1e-1,'MaxStep',.1);

%smooth observer
[t1,~,x1] = HyEQsolver(@f1,@g1,@C1,@D1,x0,TSPAN,JSPAN,rule,options);

%hybrid observer without decoupling
[t2,~,x2] = HyEQsolver(@f1,@g2,@C2,@D2,x0,TSPAN,JSPAN,rule,options);

%hybrid observer with  decoupling
[t3,~,x3] = HyEQsolver(@f2,@g2,@C2,@D2,x0,TSPAN,JSPAN,rule,options);

%hybrid observer with fully decoupling
[t4,~,x4] = HyEQsolver(@f3,@g2,@C2,@D2,x0,TSPAN,JSPAN,rule,options);
%% plot solution
Qout1     = x1(:,1:4);
Pout1     = x1(:,5:7);
Qhatout1  = x1(:,8:11);
Phatout1  = x1(:,12:14);
bwhatout1 = x1(:,15:17);
bvhatout1 = x1(:,18:20);


Qout2     = x2(:,1:4);
Pout2     = x2(:,5:7);
Qhatout2  = x2(:,8:11);
Phatout2  = x2(:,12:14);
bwhatout2 = x2(:,15:17);
bvhatout2 = x2(:,18:20);


Qout3     = x3(:,1:4);
Pout3     = x3(:,5:7);
Qhatout3  = x3(:,8:11);
Phatout3  = x3(:,12:14);
bwhatout3 = x3(:,15:17);
bvhatout3 = x3(:,18:20);


Qout4     = x4(:,1:4);
Pout4     = x4(:,5:7);
Qhatout4  = x4(:,8:11);
Phatout4  = x4(:,12:14);
bwhatout4 = x4(:,15:17);
bvhatout4 = x4(:,18:20);



angle_scale = [0 0 0 1]';
 
Thetahat1 = zeros(length(t1),1);
for i=1:length(t1)
    R     = quat2rotm(Qout1(i,:));
    Rhat1 = quat2rotm(Qhatout1(i,:));
    
 
    Thetahat1(i,:) = vrrotmat2vec(R*Rhat1')*angle_scale;   
    
end

 
Thetahat2 = zeros(length(t2),1);
for i=1:length(t2)
    R      = quat2rotm(Qout2(i,:));
    Rhat2  = quat2rotm(Qhatout2(i,:));
    
  
    Thetahat2(i)= vrrotmat2vec(R*Rhat2')*angle_scale;   
    
end

 
Thetahat3 = zeros(length(t3),1);
for i=1:length(t3)
    R      = quat2rotm(Qout3(i,:));
    Rhat3  = quat2rotm(Qhatout3(i,:));
    
 
    Thetahat3(i)= vrrotmat2vec(R*Rhat3')*angle_scale;   
    
end

 
Thetahat4 = zeros(length(t4),1);
for i=1:length(t4)
    R      = quat2rotm(Qout4(i,:));
    Rhat4  = quat2rotm(Qhatout4(i,:));
 
    Thetahat4(i)= vrrotmat2vec(R*Rhat4')*angle_scale;   
    
end


bwerror1 = sqrt(sum((bwhatout1 - ones(size(t1))*bw').^2,2));
bwerror2 = sqrt(sum((bwhatout2 - ones(size(t2))*bw').^2,2));
bwerror3 = sqrt(sum((bwhatout3 - ones(size(t3))*bw').^2,2));
bwerror4 = sqrt(sum((bwhatout4 - ones(size(t4))*bw').^2,2));

bverror1 = sqrt(sum((bvhatout1 - ones(size(t1))*bv').^2,2));
bverror2 = sqrt(sum((bvhatout2 - ones(size(t2))*bv').^2,2));
bverror3 = sqrt(sum((bvhatout3 - ones(size(t3))*bv').^2,2));
bverror4 = sqrt(sum((bvhatout4 - ones(size(t4))*bv').^2,2));


figure;
% set(Fb, 'Position', [100 100 600 560])
h1=axes('position',[0.1 0.1 0.8 0.8]);
axis(h1);
subplot(2,1,1)
h2 = plot(t1, (Thetahat1),'m-',...
     t2, (Thetahat2),'b-',...
     t3, (Thetahat3),'r-',...
     t4, (Thetahat4),'g-','linewidth',2);
ylabel('Rot.angle (rad)')
ylabel('Rot.angle (rad)','fontsize',12)
grid on
xlim([0 100])
h3=axes('position',[0.5 0.7 0.3 0.2]);
axis(h3);
h4 = plot(t1, (Thetahat1),'m-',...
     t2, (Thetahat2),'b-',...
     t3, (Thetahat3),'r-',...
     t4, (Thetahat4),'g-','linewidth',2);
xlim([0 10]),ylim([0 3.2])
grid on

subplot(2,1,2)
plot(t1,(sqrt(sum((Phatout1-Pout1).^2,2))),'m-',...
     t2,(sqrt(sum((Phatout2-Pout2).^2,2))),'b-',...
     t3,(sqrt(sum((Phatout3-Pout3).^2,2))),'r-',...
     t4,(sqrt(sum((Phatout4-Pout4).^2,2))),'g-','linewidth',2)
grid on
xlabel('time (s)')
ylabel('$\|p-\hat{p}\|$','Interpreter','latex','fontsize',14)
xlim([0 100])
legend('S','H','HD1','HD2') 


% figure
% subplot(2,1,1)
% semilogy(t1, (Thetahat1),'k-. ',...
%      t2, (Thetahat2),'b-',...
%      t3, (Thetahat3),'r-',...
%      t4, (Thetahat4),'g-','linewidth',2),
% ylabel('Rot.angle (rad)')
% ylabel('Rot.angle (rad)','fontsize',12)
% legend('Smooth observer','Non-decoupled hybrid observer','Decoupled hybrid observer','Fully decoupled hybrid observer')
% grid on
% 
% subplot(2,1,2)
% semilogy(t1,(sqrt(sum((Phatout1-Pout1).^2,2))),'k-.',...
%      t2,(sqrt(sum((Phatout2-Pout2).^2,2))),'b-',...
%      t3,(sqrt(sum((Phatout3-Pout3).^2,2))),'r-',...
%      t4,(sqrt(sum((Phatout4-Pout4).^2,2))),'g-','linewidth',2)
% grid on
% xlabel('time (s)')
% ylabel('$\|p-\hat{p}\|$','Interpreter','latex','fontsize',14)
% grid on

 


figure;
% set(Fb, 'Position', [100 100 600 560])
h1=axes('position',[0.1 0.1 0.8 0.8]);
axis(h1);
subplot(2,1,1)
h2 = plot(t1,bwerror1,'m-',t2,bwerror2,'b-',t3,bwerror3,'r-',t4,bwerror4,'g-','linewidth',2);
legend('S','H','HD1','HD2')
ylabel('$\|\tilde{b}_{\omega}\|$','Interpreter','latex','fontsize',14)
grid on
% xlim([0 100])
% % h3=axes('position',[0.45 0.3 0.43 0.35]);
% h3=axes('position',[0.45 0.7 0.4 0.2]);
% axis(h3);
% h4 = plot(t1,bwerror1,'m-',t2,bwerror2,'b-',t3,bwerror3,'r-',t4,bwerror4,'g-','linewidth',2);
% grid on

h5=axes('position',[0 0 1 1]);
axis(h5);
subplot(2,1,2)
h6 = plot(t1,bverror1,'m-',t2,bverror2,'b-',t3,bverror3,'r-',t4,bverror4,'g-','linewidth',2);
ylabel('$\|\tilde{b}_{v}\|$','Interpreter','latex','fontsize',14)
grid on
xlim([0 100])
xlabel('time (s)') 
h7=axes('position',[0.48 0.25 0.4 0.2]);
axis(h7);
h8 = plot(t1,bverror1,'m-',t2,bverror2,'b-',t3,bverror3,'r-',t4,bverror4,'g-','linewidth',2);
grid on


% figure
% plot3(Pout1(:,1),Pout1(:,2),Pout1(:,3),'k-.','linewidth',2), hold on
% plot3(Phatout1(:,1),Phatout1(:,2),Phatout1(:,3),'m-','linewidth',2), 
% plot3(Phatout2(:,1),Phatout2(:,2),Phatout2(:,3),'b-','linewidth',2),
% plot3(Phatout3(:,1),Phatout3(:,2),Phatout3(:,3),'r-','linewidth',2),
% plot3(Phatout4(:,1),Phatout4(:,2),Phatout4(:,3),'g-','linewidth',2),
% % plot3(p3(1),p3(2),p3(3),'o','MarkerSize',12),
% % plot3(p2(1),p2(2),p2(3),'o','MarkerSize',12)
% 
% xlabel('$x(m)$','Interpreter','latex','fontsize',14)
% ylabel('$y(m)$','Interpreter','latex','fontsize',14)
% zlabel('$z(m)$','Interpreter','latex','fontsize',14)
% legend('Real','Smooth-observer','Hybrid-observer','Decoupled-hybrid-observer','Fully decoupled hybrid observer')
% grid on


% figure
% subplot(2,1,1)
% plot(t1,sqrt(sum(bwhatout1.^2,2)),'m-','linewidth',2), hold on
% plot(t2,sqrt(sum(bwhatout2.^2,2)),'b-','linewidth',2)
% plot(t3,sqrt(sum(bwhatout3.^2,2)),'r-','linewidth',2)
% plot(t4,sqrt(sum(bwhatout4.^2,2)),'g-','linewidth',2)
% ylabel('$\|\hat{b}_{\omega}\|$','Interpreter','latex','fontsize',14)
% xlabel('time (s)','Interpreter','latex','fontsize',14)
% subplot(2,1,2)
% plot(t1,sqrt(sum(bvhatout1.^2,2)),'m-','linewidth',2), hold on
% plot(t2,sqrt(sum(bvhatout2.^2,2)),'b-','linewidth',2)
% plot(t3,sqrt(sum(bvhatout3.^2,2)),'r-','linewidth',2)
% plot(t4,sqrt(sum(bvhatout4.^2,2)),'g-','linewidth',2)
% ylabel('$\|\hat{b}_{v}\|$','Interpreter','latex','fontsize',14)
% xlabel('time (s)','Interpreter','latex','fontsize',14)

 

