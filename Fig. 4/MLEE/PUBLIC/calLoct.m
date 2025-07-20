function rLoct = calLoct(r, r_sw, K)
%CALLOCT Summary of this function goes here
%   Detailed explanation goes here


% r_n's location
rLoct = floor(r/r_sw)+1;
rLoct(rLoct>K) = K;

end

