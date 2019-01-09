function  xplus = g1(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file               
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% state
Q       = x(1:4,1);
p       = x(5:7,1);

Qhat    = x(8:11,1);
phat    = x(12:14,1);
bahat   = x(15:20,1);

t       = x(21,1);


 


Qhatplus = Qhat;
phatplus = phat;


xplus = [Q;p;Qhatplus;phatplus;bahat;t];




end 