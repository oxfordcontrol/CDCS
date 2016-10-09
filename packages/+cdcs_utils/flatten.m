function z = flatten(z,Z,svecFlag)

%flatten a cell array in a long vector
%For matrix cells, assume symmetric by default and use symmetric vectorization.
%To disable this behaviour and use normal matrix vectorization, set the input
%flag svecFlag to false (or 0).

if nargin>2 && svecFlag==0
    Zvec = cellfun(@(x)x(:),Z,'UniformOutput',false);
else
    Zvec = cellfun(@(x)svec(x),Z,'UniformOutput',false);
end
z = vertcat(Zvec{:});

end
