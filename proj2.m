function y=proj2(x,dx,delta)


 kb = 0;
 temp=norm(x); 

 if (temp ==0)     
     sat_delta = 0;
 else
      sat_delta = x*min(1,(delta+0.1)/temp);
 end
 
 
 y = -kb*x + kb*sat_delta + dx;
 
 

end
    
    