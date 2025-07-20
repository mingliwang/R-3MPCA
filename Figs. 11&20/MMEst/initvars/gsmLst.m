function [phi, phid] = gsmLst(gsm, N)
%GSMLST Summary of this function goes here
%   Detailed explanation goes here


% phi(s, b), s=r^2/2
idt = str2double(regexp(gsm, '(-)?\d+(\.\d+)?(e(-|+)\d+)?', 'match'));  % get numbers from char array
gsm_id = false(1, 5);
gsm_id(idt) = 1;
eps = 1e-16;
if (gsm_id(1))
    b = 5;
    phi = @(s) b*sqrt(s/N+eps);
elseif (gsm_id(2))
    b = 255^2;
    phi = @(s) log(1+b*(s/N+eps));
elseif (gsm_id(3))
    b = 20;
    phi = @(s) log(1+b*(s/N+eps));
elseif (gsm_id(4))
    b = 255;
    psi = @(s) log(1+b*sqrt(s/N+eps));
    dlt_nlf = @(s) (N-1)*log(b*sqrt(s/N+eps));
    phi = @(s) log(1+psi(s));
end


% phi(s), phid(s)
syms s
if (exist('psi', 'var'))
    lam = 1e8;
    psid = @(s) diff(psi(s));
    phi = @(s) phi(s)+(-log(psid(s))+dlt_nlf)/lam;
end
phi = matlabFunction(phi(s));
phid = matlabFunction(diff(phi(s)));

end

