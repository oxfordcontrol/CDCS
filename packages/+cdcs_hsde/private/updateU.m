function [u,others] =updateU(hatu,Y,v,others,K,rho)
%u = updateU(hatu,v,X,K,rho)
% projection to cones
    
    % cone variables: projection
    X = others.X;
    %X = cdcs_utils.blockify(X,hatu.xh-v.xh./rho,K);          % blockxified variables
    X = cdcs_utils.blockify(X,hatu.xh-v.xh,K);                % blockxified variables
    X = cdcs_utils.projectK(X,K,0);    
    u.xh = cdcs_utils.flatten(v.xh,X);
    
    % free variables
    u.x = hatu.x - v.x;%./rho;
    u.y = hatu.y - v.y;%./rho;
    u.v = hatu.v - v.v;%./rho;
    
    % non-neggative variables
    u.tau = hatu.tau - v.kappa;%./rho;
    if u.tau < 0;
        u.tau = 0;
    end

end

