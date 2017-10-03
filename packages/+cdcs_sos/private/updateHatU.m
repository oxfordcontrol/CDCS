function [hatu,others] = updateHatU(X,u,v,b,c,xi,solInner,rho,alpha,others)
% update hat{u} = (I+Q)^{-1}(u+v)

% Projection to affine set

%% v = u + v; v = u+v/rho?
% GF: rescaling v does not make any difference -> use v as the scaled multiplier
% in algorithm by O'Donoghue et at
% v.x  = u.x + (v.x)./rho;
% v.y  = u.y + (v.y)./rho;
% v.kappa  = u.tau+v.kappa./rho;
v.x  = u.x + v.x;
v.y  = u.y + v.y;
v.kappa  = u.tau + v.kappa;


%% solving (I+Q)w = v
v.x = v.x - (v.kappa).*c;
v.y = v.y + (v.kappa).*b;

w = solInner(v);

% Compute transpose only once per iteration (could use persistent variable)
ctr = c.';
btr = b.';

const    = ctr*w.x - btr*w.y;            % this is a scalar!
hatu.x   = w.x - xi.x.*const;
hatu.y   = w.y - xi.y.*const;
hatu.tau = v.kappa + ctr*hatu.x - btr*hatu.y;

%% over-relaxation; see section 3.3 in the paper by O'Donoghue
if alpha~=1
    beta = 1-alpha;
    hatu.x  = alpha.*hatu.x + beta.*u.x;
    hatu.y  = alpha.*hatu.y + beta.*u.y;
    hatu.tau = alpha.*hatu.tau + beta.*u.tau;
end


end

