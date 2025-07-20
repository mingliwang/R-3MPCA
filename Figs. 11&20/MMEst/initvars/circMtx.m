function F = circMtx(M, order, varTyp)
%CIRCMTX Summary of this function goes here
%   Detailed explanation goes here



if (M==1)
    F = 0;
    return
end


if (isempty(varTyp))
    I_M = eye(M);
else
    I_M = eye(M, varTyp);
end

order = min(max(1, order), M-1);
for k = 1: order;
    F = 0*I_M;
    for i = M:(-1):2
        F(i, :) = I_M(i, :) - I_M(i-1, :);
    end
    F(1, :) = circshift(F(2, :), [0, -1]);
    
    I_M = F;
end

F(1:order, :) = 0*F(1:order, :);
F = circshift(F, [-order, 0]);

end


