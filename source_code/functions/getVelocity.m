%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity estimation with autocorrelation method                        %
%                                                                        %
% Author: Madhavanunni A N                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function velocityMatx = getVelocity(samples_left,samples_right,aCorrLen,nAlphas,sysPara,alpha_left,bfParams)

nEns = size(samples_left,4)-aCorrLen+1;
zl = size(samples_left,2);
xl = size(samples_left,3);
nPix = xl*zl;
tx_cf = (1+bfParams.DMAS)*sysPara.fc;
% tx_cf = sysPara.fc;
constant = 2*tx_cf/sysPara.c;
for alphaIdx = 1:nAlphas
    this_samples_left = squeeze(samples_left(alphaIdx,:,:,:));
    this_samples_right = squeeze(samples_right(alphaIdx,:,:,:));
    for ensIdx = 1:nEns
        %Wall filter
        ecbf_samples_left = this_samples_left(:,:,ensIdx:ensIdx+aCorrLen-1);
        ecbf_samples_right = this_samples_right(:,:,ensIdx:ensIdx+aCorrLen-1);

%         ecbf_samples_left = wfilt(this_samples_left(:,:,ensIdx:ensIdx+aCorrLen-1),'svd',5);
%         ecbf_samples_right = wfilt(this_samples_right(:,:,ensIdx:ensIdx+aCorrLen-1),'svd',5);
        
        hilbert_ecbfLeft = reshape(hilbert(reshape(ecbf_samples_left,zl,xl*aCorrLen)),nPix,aCorrLen);
        hilbert_ecbfRight = reshape(hilbert(reshape(ecbf_samples_right,zl,xl*aCorrLen)),nPix,aCorrLen);
        
        %Autocorrelation
        aCorr_left = conj(hilbert_ecbfLeft(:,1:aCorrLen-1)).*(hilbert_ecbfLeft(:,2:aCorrLen));
        aCorr_right = conj(hilbert_ecbfRight(:,1:aCorrLen-1)).*(hilbert_ecbfRight(:,2:aCorrLen));
        
        %frequency estimation
        freqLeft = angle(mean(aCorr_left,2))*sysPara.PRF/(2*pi);
        freqRight = angle(mean(aCorr_right,2))*sysPara.PRF/(2*pi);
        f_sum = freqLeft+freqRight;
        f_sub = freqLeft-freqRight;
        
        %Velocity estimation
        cosAlpha(:,1) = cos((alpha_left(alphaIdx,:)));
        sinAlpha(:,1) = sin((alpha_left(alphaIdx,:)));
        velocityMatx.freqLeft(alphaIdx,:,ensIdx) = single(freqLeft);
        velocityMatx.freqRight(alphaIdx,:,ensIdx) = single(freqRight);
        velocityMatx.Vx(alphaIdx,:,ensIdx) = single(f_sum./((1+cosAlpha)*constant));
        velocityMatx.Vz(alphaIdx,:,ensIdx) = single(f_sub./(sinAlpha*constant));

        clc;
        disp(strcat('alphaIdx: ',num2str(alphaIdx),' --- ',num2str(ensIdx),' / ',num2str(nEns)))
    end
end

