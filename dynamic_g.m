function dx = dynamic_g(x,xi)
%% state
Q = x(1:4,1);
p = x(5:7,1);
w = xi(1:3,1);
v = xi(4:6,1);


%% differential equations
Qw = [0; w];
dQ = 0.5*quatmultiply(Q',Qw');
dQ = dQ';

R = quat2rotm(Q');
dp = R* v;

dx = [dQ;dp];