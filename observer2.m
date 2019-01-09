function dxhat = observer2(x,xi,hatba,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file              
%
% Description: Observer with dynamic decoupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global r k  k_beta
%% states
Qhat   = x(1:4,1);
phat   = x(5:7,1);
w      = xi(1:3,1);
v      = xi(4:6,1);
bawhat = hatba(1:3,1);
bavhat = hatba(4:6,1);
bahat = [bawhat;bavhat];
% Qhat = Qhat/norm(Qhat); 


%% differential equations 
[~,n] = size(r); 


Rhat = quat2rotm(Qhat');
ghat = [Rhat phat;zeros(1,3) 1]; 
ghat_inv = [Rhat' -Rhat'*phat;zeros(1,3) 1];

% find the center

pc = zeros(3,1);
temp = 0;
for i=1:n
    pc = pc + k(i)*r(1:3,i)*r(4,i);
    temp = temp + k(i)*r(4,i);
end
pc = pc/temp;
ga =  [eye(3) pc;zeros(1,3) 1];
ga_inv =  [eye(3) -pc;zeros(1,3) 1];    

rbar = zeros(size(r));
for i=1:n
   rbar(:,i) = ga_inv*r(:,i);
end
    
psi_g = zeros(6,1);
for i=1:n
%     psi_g = psi_g + 0.5*k(i)*wedge_product(ghat*b(:,i),r(:,i));
    psi_g = psi_g + 0.5*k(i)*wedge_product(ga_inv*ghat*b(:,i),rbar(:,i));
end

% k_beta = 0.5;


Ad_ghat = Ad(ghat_inv*ga);
beta = k_beta*Ad_ghat*psi_g;
beta_w = beta(1:3,1); 
beta_v = beta(4:6,1);



Ad_ghat = Ad(ga_inv*ghat);
sigma_b = Ad_ghat'*psi_g;  


 

T= [1*eye(3) zeros(3);zeros(3) 15*eye(3)];
% dbhat = proj(bahat,-T*sigma_b,0.5,T);
 
dbhat(1:3,1) = proj(bawhat,-1*sigma_b(1:3,1),0.5,1);
dbhat(4:6,1) = proj(bavhat,-10*sigma_b(4:6,1),0.5,10);


% dbhat(1:3,1) = proj(bawhat,-0.5*sigma_b(1:3,1),0.5);
% dbhat(4:6,1) = proj(bavhat,-30*sigma_b(4:6,1),0.5);



Qwhat = [0;w - bawhat + beta_w];
dQhat = 0.5*quatmultiply(Qhat',Qwhat');
dQhat = dQhat';
% dphat = Rhat*(Rhat'*v - Rhat'*bavhat + beta_v); % inertial frame velocity
dphat = Rhat*(v - bavhat + beta_v);


dxhat = [dQhat;dphat;dbhat];

 




