function adg = Ad(g)

R = g(1:3,1:3);
p = g(1:3,4);

hatp =  [0 -p(3) p(2);
         p(3) 0 -p(1);
         -p(2) p(1) 0];
 adg = [R zeros(3);hatp*R R];