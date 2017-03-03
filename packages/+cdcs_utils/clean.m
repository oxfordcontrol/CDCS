function varargout = clean(varargin)

%clean out small values in a vector.
%
% Usage : x = clean(x),     or 
%         x = clean(x,tol), or
%         clean(x,tol);

if(nargin<2)
    level  = 1e-10;
    nClean = 1;
else
    level = varargin{end};
    nClean = nargin-1;
end

%eliminate spurious values.  
%nested indexing is to avoid
%creating extremely dense
%indices from sparse matrices

%old usage
if(nargout >= 1)
    for i = 1:nClean
        varargout{i} = cleanIt(varargin{i},level);
    end
end

%super clever matlab-y amazement
if(nargout == 0)
    for i = 1:nClean
        theName = inputname(i);
        assignin('caller',theName,cleanIt(varargin{i},level));
    end
    varargout = {};
end


function x = cleanIt(x,level)

%actual cleaning function.  careful indexing
%to avoid dense indices

idx1  = find(x);
idx2  = abs(x(idx1)) < level;
x(idx1(idx2)) = 0;


