function [Pcs, Tau, Oms, nll] = Access(Y, ind_o, Rvph, Rtvs, Rtw, isTuc, X_orig, SEPs)
%ACCESS Summary of this function goes here
%   Detailed explanation goes here



% initializing
if (isTuc)
    Pcs = tucker3(Y, Rtw, true, []);
else
    Pcs = tucker1(Y, Rtw, true, []);
end


% update Pcs
LOOP = 300;
prec = 1e-3;
nll = zeros(1, LOOP);

[nll(1), Tau, Oms] = calOFV(Y, ind_o, Rvph, Rtvs, Rtw, Pcs);
for loop = 1: LOOP
    Pcs = P_step(Y, Rvph, Rtvs, Rtw, Pcs, Tau, Oms);
    
    [nll(loop+1), Tau, Oms] = calOFV(Y, ind_o, Rvph, Rtvs, Rtw, Pcs);
    if (abs(nll(loop+1)-nll(loop)) <= prec)
        break
    end
    
    imDisp(Y, ind_o, Rtw, Pcs, loop, X_orig, SEPs);
    
end

nll = nll(1:(1+loop));                                            % target function value

end
