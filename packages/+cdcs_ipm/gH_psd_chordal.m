function [in, g, H, L] = gH_psd_chordal(x, params)
    %sum of all the log-det penalties in each clique
    %I don't know if this is more or less efficient than the Vandenberghe
    %scheme of recursing up and down the clique tree, especially if the
    %penalty of their scheme is only applied to one side.
    %x_clique = x(params.Ech);
    params_inner = params;
    NoElem = params.NoElem;
    
    g = zeros(size(x));
    H = sparse(length(x), length(x));
    
    %g_accum = 
    %H_accum = zeros(size(params.Ech));
    %for each clique, apply log det penalty
    count = 0;
    in = 1;
    %H_cell = {}
    for i = 1:length(NoElem)        
        %find contents of each clique
        NoVars = NoElem(i)*(NoElem(i) + 1)/2;
        ind = count+1:count+NoVars;
        ind_Ech = params.Ech(ind);
        x_k = x(ind_Ech); 
        
        %get gradient and hessian
        params_inner.Q = params.Q(NoElem(i));
        [in_k, g_k, H_k] = gH_psd_inner(x_k, params_inner);
        in = in && in_k;
        
        %accumulate them
        g(ind_Ech) = g(ind_Ech) + g_k;
        H(ind_Ech, ind_Ech) = H(ind_Ech, ind_Ech) + H_k;
        %H_accum(ind) = H_k;
        count = count + NoVars;
        %H_cell{i}.Ech = ind_Ech;
        %H_cell{i}.H   = H_k;
    end
    
    %H_raw = accumarray(Ech,H_accum, [], sum);

    %L = chol(H, 'lower');
    if in
        [L, err] = chol(H, 'lower');
        %these will be sparse cholesky factors, take advantage of this.
        %[L, err, R] = chol(H, 'lower');
        in  = in & (err == 0);
    else
        g = NaN;
        H = NaN;
        L = NaN;
    end
end