function y=proj(x,dx,delta,T)


e = 0.1;
Px = norm(x) - delta;

temp = norm(x);
if (temp ==0)
    dPx=zeros(size(x));
else
    dPx = x/(temp);
end
% Px = 0.5*(x'*x - delta^2);
% dPx = x;

c = max(1,Px/e);
% c=1;

y = zeros(size(dx));

% if (Px <=0)||((dPx'*dx < 0)&&((Px<=e)&&(0<Px)))
if (Px <= 0)||(dPx'*dx <=0)
    y = dx;
else 
% elseif ((dPx'*dx > 0)&&((Px<=e)&&(0<Px)))
%     y = (eye(length(x))-(x*x'/(x'*x)))*dx;
%     Px

    y = (eye(length(x))-c*T*(dPx*dPx'/(dPx'*T*dPx)))*dx;
    

end
% norm((eye(length(x))-c*(dPx*dPx'/(dPx'*dPx))))


end
    
    