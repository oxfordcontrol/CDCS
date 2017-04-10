function [v_n,others]=updateV(hatu,u,v,rho,others)
%   v_n = updateV(v,hatu,u,rho)
%   The third step in ADMM, update the multipliers

%     v_n.x  = v.x - rho.*(hatu.x - u.x);       %rho*(hatu.x - u.x);
%     v_n.y  = v.y - rho.*(hatu.y - u.y);       %rho*(hatu.y - u.y);
%     v_n.kappa = v.kappa - rho.*(hatu.tau - u.tau);%rho*(hatu.kappa - u.kappa);

% GF: v is already scaled multiplier
v_n.x  = v.x - hatu.x + u.x;       
v_n.y  = v.y - hatu.y + u.y;      
v_n.kappa = v.kappa - hatu.tau + u.tau;


end

