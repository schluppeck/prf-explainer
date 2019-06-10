function [K] = estimate_prf_linear_transform(stimImage, hrf)
%estimate_prf_linear_transform - estimate RF via linear regression methods
%   
% .   this is a helper function that sets up the linear problem
% .   and returns the matrix K that needs to be inverted with
% .   regularisation
% .   
%
% ds 2019-06-04


% S is the stimulus (
% .       rows: space unwrapped into a numel(pRF) vector)
% .       columns: time points
% .       )
% H is the convolution matrix (HRF) 

% K = H * S

% making H:

nTimepoints = size(stimImage.im,3);
H = convmtx(hrf, nTimepoints);
H = H(1:nTimepoints, 1:nTimepoints);

% making S
% first two dimenions of stimImage.im: image size;, third: time
nPixels = numel(stimImage.im(:,:,1));   

% "vectorise" stimImage and also the pRF
S = reshape(stimImage.im, nPixels, []);

% the K matrix is the combination of stimulus and haemodynamics
K = H * S';


end
    
    