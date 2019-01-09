function  xplus = g2(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file               
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global  AA index U theta_q

% state
Q       = x(1:4,1);
p       = x(5:7,1);

Qhat    = x(8:11,1);
phat    = x(12:14,1);
bahat   = x(15:20,1);

t       = x(21,1);


B = AA(1:3,4);
d = AA(4,4);
 
if (index ==0)
    Qq_inv = [1 0 0 0]';
    pq = 0;
else
    uq = U(:,index);
    Qq = [cos(theta_q/2) sin(theta_q/2)*uq']';
    Rq = expm(theta_q*skew(uq));
    pq = (eye(3)-Rq)*(1/d)*B;
    Qq_inv = quatinv(Qq');
    Qq_inv = Qq_inv';
end


Qhatplus = quatmultiply(Qq_inv',Qhat');
Qhatplus = Qhatplus';
phatplus = Rq'*(phat - pq);

xplus = [Q;p;Qhatplus;phatplus;bahat;t];




end 