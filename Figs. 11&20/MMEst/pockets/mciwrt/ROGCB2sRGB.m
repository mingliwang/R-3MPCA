function sRGB = ROGCB2sRGB(ROGCB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date     : Nov. 12, 2015
% Author   : Yusuke Monno
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ROGCB = clip(ROGCB,0,1);

R = ROGCB(:, :, 1);
O = ROGCB(:, :, 2);
G = ROGCB(:, :, 3);
C = ROGCB(:, :, 4);
B = ROGCB(:, :, 5);


% degamma correction
gam = 2.2;
R = degamma_correction(R,gam);
O = degamma_correction(O,gam);
G = degamma_correction(G,gam);
C = degamma_correction(C,gam);
B = degamma_correction(B,gam);



% Image size
[h,w] = size(R);

% Transformation matrix
M = [1.342586196, 3.303167555, -1.002982634, -0.49354491, 0.74854673;
    -0.043234884, -0.880120819, 1.78350697,	-0.314556273, -0.24018848;
    -0.17595671, 0.216464157, -0.274041435, -0.211740902, 2.878603932];

% Vector form of input 5-band images
vec = [R(:),O(:),G(:),C(:),B(:)]';

% Transformation
sRGB = M * vec;

% Reshpae to image form
sRGB = reshape(sRGB',h,w,3);

% Clip to [0,1]
sRGB = clip(sRGB,0,1);



% gamma correction
sRGB = gamma_correction(sRGB,gam);


% to [0, 255]/255
sRGB = uint8(255*sRGB);
sRGB = cast(sRGB, 'like', ROGCB)/255;

end


function Y = degamma_correction(X,gam)
%
%

Y = X.^gam;

end


function Y = gamma_correction(X,gam)
%
%

Y = X.^(1/gam);

end


function [Y] = clip(X,lo,hi)
%
%

Y = X;
Y(X<lo)=lo;
Y(X>hi)=hi;

end


