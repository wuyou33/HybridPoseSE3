function [value] = C1(x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file                
%
% Description: Flow set
% Return 0 if outside of C, and 1 if inside C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global  AA DD U theta_q delta

Q       = x(1:4,1);
p       = x(5:7,1);
Qhat    = x(8:11,1);
phat    = x(12:14,1);
t       = x(21,1);

% find the center
% pc = zeros(3,1);
% temp = 0;
% [~,n] = size(r);
% for i=1:n
%     pc = pc + k(i)*r(1:3,i)*r(4,i);
%     temp = temp + k(i)*r(4,i);
% end
% pc = pc/temp;


B = AA(1:3,4);
d = AA(4,4);

ga = [eye(3) B/d;zeros(1,3) 1];
ga_inv = [eye(3) -B/d;zeros(1,3) 1];


R  = quat2rotm(Q');
Rhat= quat2rotm(Qhat');
g = [R,p; zeros(1,3) 1];
ghat_inv = [Rhat', -Rhat'*phat; zeros(1,3) 1];


gtilde = ga_inv*g*ghat_inv*ga;





U_0 = 0.5*trace((eye(4)-gtilde)*DD*(eye(4)-gtilde)');


u1 = U(:,1);
R1 = expm(theta_q*skew(u1));
p1 = zeros(3,1);
g1 = [R1 p1;zeros(1,3) 1];

U_q(1)= 0.5*trace((eye(4)-gtilde*g1)*DD*(eye(4)-gtilde*g1)');

u2 = U(:,2);
R2 = expm(theta_q*skew(u2));
p2 = zeros(3,1);
g2 = [R2 p2;zeros(1,3) 1];

U_q(2) = 0.5*trace((eye(4)-gtilde*g2)*DD*(eye(4)-gtilde*g2)');

u3 = U(:,3);
R3 = expm(theta_q*skew(u3));
p3 = zeros(3,1);
g3 = [R3 p3;zeros(1,3) 1];

U_q(3) = 0.5*trace((eye(4)-gtilde*g3)*DD*(eye(4)-gtilde*g3)');

 
 

min_U = min(U_q);

if (U_0-min_U <= delta)
    value = 1;
else 
    value = 0;
end

% value = 1;

end