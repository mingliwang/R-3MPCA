function Jm = imNoise(Im, arg3, noisTyp)
%IMNOISE Noise for color video.
%

Im = single(Im);
[R, C, D, P] = size(Im);

switch noisTyp
    case 'Gaussian'
        randn('seed', 0);  % rng(0, 'v4') gives same result without warning
        noise = randn([R, C, D, P]);
        noise = permute(noise, [1, 2, 4, 3]);
        noise = reshape(noise, [], D)*diag(arg3);
        noise = reshape(noise, [R, C, P, D]);
        Jm = Im + permute(noise, [1, 2, 4, 3]);
end
end

