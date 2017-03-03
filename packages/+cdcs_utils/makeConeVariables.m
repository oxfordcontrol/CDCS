function [varargout] = makeConeVariables(K)

% Make conic variables
varargout = repmat({{}},[1 nargout]);
shift = 0;

%populate the free variables
if(isfield(K,'f') && K.f > 0)
    for j = 1:nargout
        varargout{j}{1+shift} = zeros(K.f,1);
    end
    shift = 1;
end

%populate the zero variables
if(isfield(K,'z') && K.z > 0)
    for j = 1:nargout
        varargout{j}{1+shift} = zeros(K.z,1);
    end
    shift = shift+1;
end

%populate the linear cone variables
if(isfield(K,'l') && K.l > 0)
    for j = 1:nargout
        varargout{j}{1+shift} = zeros(K.l,1);
    end
    shift = shift+1;
end

%populate the second order cone variables
if isfield(K,'q') && sum(K.q) > 0
    nqCones = length(K.q);
    for i = 1:nqCones
        for j = 1:nargout
            varargout{j}{end+1} = zeros(K.q(i),1);
        end
    end
end

%populate the semidefinite cones
if isfield(K,'s') && sum(K.s) > 0
    nsCones = length(K.s);
    for i = 1:nsCones
        for j = 1:nargout
            varargout{j}{end+1} = zeros(K.s(i),K.s(i));
        end
    end
end
