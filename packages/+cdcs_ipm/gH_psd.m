function [in, g, H, L] = gH_psd(x, params)
    %gradient and hessian of log-det penalty (with cholesky factorization)
    
    [in, g, H] = gH_psd_inner(x, params);
    
    [L, err] = chol(H, 'lower');
    in  = in & (err == 0);
    
end