function X = tnsMult(S, U_1, U_2, V)
%TNSMULT Summary of this function goes here
%   Detailed explanation goes here



% size, absent_mde
SIZ = tnsSize(S);
abs_1 = isempty(U_1);
abs_2 = isempty(U_2);
abs_3 = isempty(V);



% mode-1, 2
Z = S;                                                                 % no mode-1, 2
if (abs_2)
    if (~abs_1)
        Z = reshape(S, SIZ(1), SIZ(2)*SIZ(3));
        SIZ(1) = size(U_1, 1);
        Z = reshape(U_1*Z, SIZ);                           % only mode-1
    end
    
else
    U_2t = U_2';
    if (abs_1)
        SIZ(2) = size(U_2, 1);
        Z = zeros(SIZ, 'like', S);
        for k = 1: SIZ(3)
            Z(:, :, k) = S(:, :, k)*U_2t;                        % only mode-2
        end
    else
        SIZ(1) = size(U_1, 1);
        SIZ(2) = size(U_2, 1);
        Z = zeros(SIZ, 'like', S);
        for k = 1: SIZ(3)
            Z(:, :, k) = U_1*S(:, :, k)*U_2t;                % mode-1, 2
        end
    end
end



% mode-3
X = Z;
if (~abs_3)
    Z = reshape(Z, SIZ(1)*SIZ(2), SIZ(3));
    SIZ(3) = size(V, 1);
    X = reshape(Z*(V'), SIZ);
end

end

