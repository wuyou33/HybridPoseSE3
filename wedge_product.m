function x1x2 = wedge_product(x1,x2)


x1_v = x1(1:3,1);
x1_s = x1(4,1);

x2_v = x2(1:3,1);
x2_s = x2(4,1);

x1x2(1:3,1) = cross(x1_v,x2_v);
x1x2(4:6,1) = x1_s * x2_v - x2_s * x1_v;