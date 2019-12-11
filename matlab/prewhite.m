function [noiseMcPreWhite] = prewhite(noiseMc,apprxNoisePsd,apprxCleanPsd,mcSz,patchSz,R)
% prewhite the micrograph using the noise RPSD
%  
%  Input:
%   noiseMc            the micrograph we want to prewhiten.
%   apprxNoisePsd      RPSD of the noise
%   apprxCleanPsd      RPSD of the particle
%   patchSz            The patch size from whom we estimate the noise RPSD.
%   R                  The radii samples of the RPSD
%  Output:
%   noiseMcPreWhite    the prewhitened micrograph           
% 
%
% Amitay Eldar November 2019.

L = floor((mcSz-1)/2);
T = 1; % Nyquist sampling rate
bandLimit = pi/T;
x = (-L:1:L)*(bandLimit/L);
[X,Y] = meshgrid(x);
radMat = sqrt(X.^2+Y.^2);
[radSamp,~,idx] = unique(radMat);
noisePsdNodes = abs(spline(R*bandLimit,apprxNoisePsd,radSamp));
noisePsdMat = reshape(noisePsdNodes(idx),length(x),length(x)); % This mat will use for whitening.
noiseMcPreWhite = cryo_prewhiten(noiseMc , noisePsdMat); % we dont want zeros for whitening apprxCleanPsd.
% apprxCleanPsd = (apprxCleanPsd./apprxNoisePsd)*(sum(apprxNoisePsd)); % clean psd after whitening
% apprx new noiseVar
% stedMatGpu = stdfilt(gpuArray(noiseMc),ones(patchSz));
% stedMat = gather(stedMatGpu);
% varMat = stedMat.^2;
% cut = floor((patchSz-1)/2) + 1; % we take in considiration only whole blocks
% varMat = varMat(cut:end-cut,cut:end-cut);
% varVec = sort(varMat(:));
% j = floor(0.25*length(varVec));
% noiseVar = mean(varVec(1:j));