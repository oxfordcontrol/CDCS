function hatu = updateHatU(u,v,b,c,xi,solInner,rho,alpha)
% update hat{u} = (I+Q)^{-1}(u+v)

% GF: when multiplying a vector v by a scalars , using s.*v is usually quicker
% than s*v. Also, if v is a real vector (as I assume is the case here), it is
% better to transpose with v.' instead of v' (avoid conjugation).

% Projection to affine set

    %rho = opts.rho;
    %alpha % over relaxation?
    
    %% v = u + v; v = u+v/rho?
    v.x  = u.x + (v.x)./rho;
    v.xh = u.xh +(v.xh)./rho;
    v.y  = u.y + (v.y)./rho;
    v.v  = u.v + (v.v)./rho;
    v.kappa  = u.tau+v.kappa./rho;
    
    %% solving (I+Q)w = v
    v.x = v.x - (v.kappa).*c;
    v.y = v.y + (v.kappa).*b;
    
    w = solInner(v);
    
    % GF:
    % Compute transpose only once per iteration (could use persistent variable)
    ctr = c.';
    btr = b.';
    
    const    = ctr*w.x - btr*w.y;            % this is a scalar!
    hatu.x   = w.x - xi.x.*const;
    hatu.xh  = w.xh - xi.xh.*const;
    hatu.y   = w.y - xi.y.*const;
    hatu.v   = w.v - xi.v.*const;  
    hatu.tau = v.kappa + ctr*hatu.x - btr*hatu.y;
    
    %% over-relaxation; see section 3.3 in the paper by O'Donoghue
    if alpha~=1
        % GF:
        % compute only once per iteration (in fact, could compute only once and
        % for all if alpha does not change between iterations!)
        beta = 1-alpha;                        
        hatu.x  = alpha.*hatu.x + beta.*u.x;
        hatu.xh = alpha.*hatu.xh + beta.*u.xh;
        hatu.y  = alpha.*hatu.y + beta.*u.y;
        hatu.v  = alpha.*hatu.v + beta.*u.v;  
        hatu.tau = alpha.*hatu.tau + beta.*u.tau;
    end


end

