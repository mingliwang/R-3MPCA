function Rvph = initVph(r_dp, gsm, isBTF, N)
%INITVPH Summary of this function goes here
%   Detailed explanation goes here



if (isBTF)
    r_dp = r_dp*sqrt(N);
else
    N = cast(1, 'like', N);
end
[phi, phid] = gsmLst(gsm, N);


% h(r) = tau_1*(r^2/2)*I_1 + (phi(r^2/2)+beta)*I_2
% I_1=(0<r<=r_dp), I_2=I_2(r>r_dp)
s_dp = (r_dp^2)/2;
tau_1 = phid(s_dp);
beta_phi = tau_1*(s_dp) - phi(s_dp);


if (isBTF)
    Rvph.isBTF = isBTF;
    Rvph.replc = 0.1;                                           % replaced cons, for missing element
end
Rvph.phi = phi;
Rvph.phid = phid;
Rvph.r_dp = r_dp;
Rvph.tau_1 = tau_1;
Rvph.beta_phi = beta_phi;


Rvph.theta = cast(1, 'like', r_dp);                        % fixed; Rpsi.mu, Rwv.regs

end

