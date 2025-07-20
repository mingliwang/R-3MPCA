function [Y_etr, ind_otr, SEPs] = tnsStre(Y, ind_o, scals)
%TNSSTRE Summary of this function goes here
%   Detailed explanation goes here




SIZ = tnsSize(Y);
M_1 = SIZ(1);
M_2 = SIZ(2);

SCL = 2^scals;
if (mod(SIZ(1), SCL)~=0)
    M_1 = (floor(SIZ(1)/SCL)+1)*SCL;
end
if (mod(SIZ(2), SCL)~=0)
    M_2 = (floor(SIZ(2)/SCL)+1)*SCL;
end


[Y_etr, SEPs] = stretch(Y, SIZ, M_1, M_2);
[ind_otr, ~] = stretch(ind_o, SIZ, M_1, M_2);

end



function [Y_etr, SEPs] = stretch(Y, SIZ, M_1, M_2)
%STRETCH Summary of this function goes here
%   Detailed explanation goes here



part_r1 = [];
part_r2 = [];
if (SIZ(1)~=M_1)
    R_1 = floor((M_1-SIZ(1))/2);
    R_2 = M_1-SIZ(1)-R_1;
    
    part_r1 = Y(1:R_1, :, :);
    part_r2 = Y((SIZ(1)-R_2+1):SIZ(1), :, :);
    
    part_r1 = flip(part_r1, 1);
    part_r2 = flip(part_r2, 1);
    
else
    R_1 = 0;
    
end
Y_etr = cat(1, part_r1, Y, part_r2);


part_c1 = [];
part_c2 = [];
if (SIZ(2)~=M_2)
    C_1 = floor((M_2-SIZ(2))/2);
    C_2 = M_2-SIZ(2)-C_1;
    
    part_c1 = Y_etr(:, 1:C_1, :);
    part_c2 = Y_etr(:, (SIZ(2)-C_2+1):SIZ(2), :);
    
    part_c1 = flip(part_c1, 2);
    part_c2 = flip(part_c2, 2);
    
else
    C_1 = 0;
    
end
Y_etr = cat(2, part_c1, Y_etr, part_c2);


SEPs(1, 1, 1) = R_1+1;
SEPs(1, 2, 1) = R_1+SIZ(1);
SEPs(1, 1, 2) = C_1+1;
SEPs(1, 2, 2) = C_1+SIZ(2);

end



