function v_n = updateV(v,hatu,u,rho)
%   The third step in ADMM, update the multipliers

    v_n.x  = v.x - rho.*(hatu.x - u.x);       %rho*(hatu.x - u.x);
    v_n.xh = v.xh - rho.*(hatu.xh - u.xh);    %rho*(hatu.xh - u.xh);
    v_n.y  = v.y - rho.*(hatu.y - u.y);       %rho*(hatu.y - u.y);
    v_n.v = v.v - rho.*(hatu.v - u.v);        %rho*(hatu.v - u.v);
    v_n.kappa = v.kappa - rho.*(hatu.tau - u.tau);%rho*(hatu.kappa - u.kappa);
    
end

