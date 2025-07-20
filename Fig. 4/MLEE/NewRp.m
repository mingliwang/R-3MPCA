function varsigj = NewRp(r_e, M, K, A_l, A_u, h_grad, Z_pt1, Z_pp1, prec, Z_ign, varsig, j)
%NEWRP Newton-Raphson method for updating varsigj
%   Detailed explanation goes here.



% tau, rho
tau = A_u*varsig;
rho = (A_l*varsig)';

% Z_pt
Z_pt = zeros(3, K, 'double');
Z_pt(1, :) = Z_pt1;  % Z_pt(1, :)=-(A_l(:, j)+(1/2)*A_u(:, j))'; k=1,...,j, then k=j+1,...,K;
Z_pt(2, :) = exp(-rho);
M_pt = M + 2*A_u(:, j);
for k = 1: K
    Z_pt(3, k) = gmsinc(r_e(k), r_e(k+1), M_pt(k), tau(k));
end

% Z_pp
Z_pp = Z_pt;
Z_pp(1, :) = Z_pp1;  % Z_pp(1, :) = ((A_l(:, j)).^2+(1/4)*A_u(:, j))';
M_pp = M + 4;
for k = 1: j
    Z_pp(3, k) = gmsinc(r_e(k), r_e(k+1), M_pp, tau(k));
end


% Newton-Raphson method
CIRC = 1e3;
varsigj = varsig(j);
for circ = 1: CIRC
    % Z_de, Z_dd
    Z_de = sum((Z_pt(1, :).*Z_pt(2, :)).*Z_pt(3, :));
    Z_dd = sum((Z_pp(1, :).*Z_pp(2, :)).*Z_pp(3, :));
    
    % G_de, G_dd
    G_de = Z_ign*Z_de + h_grad(j);
    G_dd = Z_ign*Z_dd;
    
    
    % varsigj_upd
    decre = G_de/G_dd;
    if (varsigj<=decre)
        decre = varsigj/2;
    end
    varsigj_upd = varsigj - decre;
    
    if (abs(decre)<prec)
        break
    else
        varsigj = varsigj_upd;
        tau(1:j) = tau(1:j) - decre;  % tauj=varsigj, tauj_upd=varsigj_upd;
        rho((j+1):K) = rho((j+1):K) - A_l(K, j)*decre;  % rho=A_l*(varsigj-decre);
    end
    
    
    % Z_pt, Z_pp
    Z_pt(2, (j+1):K) = exp(-rho((j+1):K));
    Z_pp(2, (j+1):K) = Z_pt(2, (j+1):K);
    for k = 1: j
        Z_pt(3, k) = gmsinc(r_e(k), r_e(k+1), M_pt(k), tau(k));
        Z_pp(3, k) = gmsinc(r_e(k), r_e(k+1), M_pp, tau(k));
    end
    
end
end