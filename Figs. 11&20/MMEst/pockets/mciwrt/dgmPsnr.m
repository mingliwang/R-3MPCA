function p = dgmPsnr(x, y)
%DGMPSNR See Sensors2017_ARI_Code
%


% [0, 255]/255
x = double(uint8(255*x))/255;
y = double(uint8(255*y))/255;


% degamma correction
gam = 2.2;
x = degamma_correction(x, gam);
y = degamma_correction(y, gam);


% psnr
[M, N, C] = size(x);
p = zeros(1, C);

err = x - y;
for c = 1: C
    e = err(:, :, c);
    mse=sum(sum(e.*e))/(M*N);
    
    if (mse>0)
        p_c=10*log(1/mse)/log(10);
    else
        p_c=99.99;
    end
    p(c) = p_c;
end

p = mean(p(:));

end



function Y = degamma_correction(X,gam)
%
%

Y = X.^gam;

end



