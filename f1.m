function xdot = f1(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file              
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global r ba delta_m

% state
Q       = x(1:4,1);
p       = x(5:7,1);
Qhat    = x(8:11,1);
phat    = x(12:14,1);
bawhat  = x(15:17,1);
bavhat  = x(18:20,1);
t       = x(21,1);

%% measurements
% generate velocity input
w    = 1*[-sin(20*t) cos(20*t) 0]';
v    = 15*[cos(0.1*t) sin(0.1*t) 0]';
xi   = [w;v];

bavary = ba;%*sin(0.001*t);
xi_y = [w;v] + bavary  ;  %baised velocity

R = quat2rotm(Q');
g_inv = [R' -R'*p; zeros(1,3) 1];
[~,n] = size(r);
b=zeros(size(r));
for i=1:n
   b(:,i) = g_inv*r(:,i) + delta_m*[randn(3,1);0];
end

%% differential equations

x_now = [Q;p];
% call real pose dynamics
dx = dynamic_g(x_now,xi);

% xi_y = [w;R*v] + ba  ;  %baised velocity
% call observer
hatx = [Qhat;phat];
hatba  = [bawhat;bavhat];
dxhat = observer1(hatx,xi_y,hatba,b);


xdot = [dx;dxhat;1];

end