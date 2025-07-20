function c = gmsinc(a, b, M, tau)
%GMSINC Summary of this function goes here
%   Integral from a to b of r^(M-1) *exp(-tau*(r^2)/2) dr.



a = double((tau/2)*(a^2));
b = double((tau/2)*(b^2));

c_a = gammainc(a, M/2, 'upper');
c_b = gammainc(b, M/2, 'upper');

c = ((2/tau)^(M/2))*(gamma(M/2)/2)*(abs(c_b-c_a));
end

