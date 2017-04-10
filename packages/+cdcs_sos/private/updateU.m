function [u,others] =updateU(hatu,Y,v,others,K,rho)
%u = updateU(hatu,v,X,K,rho)
% projection to cones

% cone variables: projection
X = others.X;
% GF: make v scaled multiplier
% X = cdcs_utils.blockify(X,hatu.x-v.x./rho,K);            
X = cdcs_utils.blockify(X,hatu.x-v.x,K);              % blockified variables
X = cdcs_utils.projectK(X,K,0);
u.x = cdcs_utils.flatten(v.x,X);

% free variables
% u.y = hatu.y - v.y./rho;
u.y = hatu.y - v.y;

% non-neggative variables
% u.tau = hatu.tau - v.kappa./rho;
u.tau = hatu.tau - v.kappa;
if u.tau < 0
    u.tau = 0;
end

end

