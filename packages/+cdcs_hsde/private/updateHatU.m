function [hatu,others] = updateHatU(X,u,v,b,c,btr,ctr,xi,solInner,rho,alpha,others)
% update hat{u} = (I+Q)^{-1}(u+v)

% Projection to affine set
   
    %% v = u + v; v = u+v/rho?
%     v.x  = u.x + (v.x)./rho;
%     v.xh = u.xh +(v.xh)./rho;
%     v.y  = u.y + (v.y)./rho;
%     v.v  = u.v + (v.v)./rho;
%     v.kappa  = u.tau+v.kappa./rho;
    
    v.x  = u.x  + v.x;
    v.xh = u.xh + v.xh;
    v.y  = u.y  + v.y;
    v.v  = u.v  + v.v;
    v.kappa  = u.tau+v.kappa;
    
    %% solving (I+Q)w = v
    v.x = v.x - (v.kappa).*c;
    v.y = v.y + (v.kappa).*b;
    
    w = solInner(v);
    
    const    = ctr*w.x - btr*w.y;            % this is a scalar!
    hatu.x   = w.x - xi.x.*const;
    hatu.xh  = w.xh - xi.xh.*const;
    hatu.y   = w.y - xi.y.*const;
    hatu.v   = w.v - xi.v.*const;  
    hatu.tau = v.kappa + ctr*hatu.x - btr*hatu.y;
    
    %% over-relaxation; see section 3.3 in the paper by O'Donoghue
    if alpha~=1     
        beta = 1-alpha;  
        hatu.x  = alpha.*hatu.x + beta.*u.x;
        hatu.xh = alpha.*hatu.xh + beta.*u.xh;
        hatu.y  = alpha.*hatu.y + beta.*u.y;
        hatu.v  = alpha.*hatu.v + beta.*u.v;  
        hatu.tau = alpha.*hatu.tau + beta.*u.tau;
    end


end

