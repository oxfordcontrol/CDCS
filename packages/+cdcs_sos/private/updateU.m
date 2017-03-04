function [u,others] =updateU(hatu,Y,v,others,K,rho)
%u = updateU(hatu,v,X,K,rho)
% projection to cones
    
    % cone variables: projection
    X = others.X;
    X = cdcs_utils.blockify(X,hatu.x-v.x./rho,K);              % blockxified variables
    X = cdcs_utils.projectK(X,K,0);    
    u.x = cdcs_utils.flatten(v.x,X);
    
    % free variables
    u.y = hatu.y - v.y./rho;
    
    % non-neggative variables
    u.tau = hatu.tau - v.kappa./rho;
    if u.tau < 0;
        u.tau = 0;
    end

end

