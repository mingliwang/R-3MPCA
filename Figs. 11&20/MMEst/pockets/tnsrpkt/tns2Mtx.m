function X_mde = tns2Mtx(X, idtMde)
%TNS2MTX Tensor to matrix
%   Detailed explanation goes here



[M_1, M_2, N] = size(X);

if (idtMde(1))
    X_mde = reshape(X, M_1, M_2*N);
    
elseif (idtMde(2))
    X_mde = permute(X, [2, 1, 3]);
    X_mde = reshape(X_mde, M_2, M_1*N);
    
elseif (idtMde(3))
    X_mde = reshape(X, M_1*M_2, N);
    X_mde = X_mde';
    
end

end

