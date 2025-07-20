function X = mtx2Tns(X_mde, SIZ, idtMde)
%MTX2TNS Matrix to tensor
%   Detailed explanation goes here



if ((prod(SIZ)~=numel(X_mde)) || (length(SIZ)~=3))
    return
end


X = zeros(SIZ, 'like', X_mde);
if (idtMde(1))
    X = reshape(X_mde, SIZ);
    
elseif (idtMde(2))
    X = reshape(X_mde, SIZ(2), SIZ(1), SIZ(3));
    X = permute(X, [2, 1, 3]);
    
elseif (idtMde(3))
    X_mde = X_mde';
    X = reshape(X_mde, SIZ);
    
end

end

