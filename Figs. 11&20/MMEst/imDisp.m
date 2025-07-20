function [] = imDisp(Y, ind_o, Rtw, Pcs, loop, X_orig, SEPs)
%IMDISP Summary of this function goes here
%   Detailed explanation goes here



if (isempty(SEPs))
    X_hat = tnsMult(Pcs.S, Pcs.U_1, Pcs.U_2, Pcs.V);
else
    X_hat = tnsMult(wvSyns(Pcs.S, Rtw), [], [], Pcs.V);
    X_hat = tnsTrm(X_hat, SEPs);
    
    Y = tnsTrm(Y, SEPs);
    ind_o = tnsTrm(ind_o, SEPs);
end



if (~iscell(X_orig))
    if (size(X_orig, 3)==1)
        mae_obs = 255*mean(abs(Y(:)-X_orig(:)));
        mae_val = 255*mean(abs(X_hat(:)-X_orig(:)));
        X_hat1 = X_hat;
    else
        mae_obs = 255*mean(mean(abs(Y(:, :, 5)-X_orig(:, :, 5))));
        mae_val = 255*mean(mean(abs(X_hat(:, :, 5)-X_orig(:, :, 5))));
        X_hat1 = X_hat(:, :, 5);
        Y = Y(:, :, 5);
        ind_o = ind_o(:, :, 5);
    end
    
    pause(1e-5)
    subplot(1,2,1);
    imshow(uint8(255*(Y.*ind_o)));
    set(gca, 'Position', [0.1, 0.2, 0.4, 0.6]);
    text(0.5, -0.1, sprintf('sparse noise level 0.6 (genNm.m) \n mae = %.4f', mae_obs), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    
    subplot(1,2,2);
    imshow(uint8(255*X_hat1));
    set(gca, 'Position', [0.55, 0.2, 0.4, 0.6]);
    if (strcmpi(class(Y), 'single'))
        if (size(X_orig, 3)==10)
            text(0.5, -0.1, sprintf('Ours_2 (ten images as input): \n iteration = %d, mae = %.4f', loop, mae_val), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
        elseif (size(X_orig, 3)==1)
            text(0.5, -0.1, sprintf('Ours_2 (single image as input): \n iteration = %d, mae = %.4f', loop, mae_val), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
        end
    else
        text(0.5, -0.1, sprintf('Ours_3 (based on wavelet): \n iteration = %d, mae = %.4f', loop, mae_val), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    
else
    X_hat3 = ROGCB2sRGB(X_hat);
    psnr3 = psnr(uint8(255*X_hat3), uint8(255*X_orig{2}));
    
    pause(1e-5)
    subplot(1,2,1);
    Y = ROGCB2sRGB(Y);
    psnr3y = psnr(uint8(255*Y), uint8(255*X_orig{2}));
    imshow(uint8(255*Y));
    set(gca, 'Position', [0.1, 0.2, 0.4, 0.6]);
    text(0.5, -0.1, sprintf('sparse noise level 0.6 (genNm.m) \n psnr = %.4f', psnr3y), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    
    subplot(1,2,2);
    imshow(uint8(255*X_hat3));
    set(gca, 'Position', [0.55, 0.2, 0.4, 0.6]);
    text(0.5, -0.1, sprintf('Ours_2: \n iteration = %d, psnr = %.4f', loop, psnr3), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
end

end

