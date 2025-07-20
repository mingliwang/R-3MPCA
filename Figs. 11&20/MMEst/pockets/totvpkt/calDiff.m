function G = calDiff(X, SIZ, F, idtMde)
%CALDIFF Summary of this function goes here
%   Detailed explanation goes here



G = 0*X;
if (idtMde(1))
    X = reshape(X, SIZ(1), prod(SIZ(2:3)));
    G = reshape(F*X, SIZ);
    
elseif (idtMde(2))
    Ft = F';
    for n = 1: SIZ(3)
        G(:, :, n) = X(:, :, n)*Ft;
    end
    
elseif (idtMde(3))
    X = reshape(X, prod(SIZ(1:2)), SIZ(3));
    G = reshape(X*(F'), SIZ);
    
end

end
