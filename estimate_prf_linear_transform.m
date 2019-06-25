function [K, H] = estimate_prf_linear_transform(stimImage, hrf)
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
H = convmtx(hrf, nTimepoints+numel(hrf)-1);
offset = numel(hrf);
idx = offset + (1:nTimepoints) - 1;
H = transpose( H(idx, idx) ); % for consistency with equations (conv in time!)

% making S
% first two dimenions of stimImage.im: image size;, third: time
nPixels = numel(stimImage.im(:,:,1));   

% "vectorise" stimImage and also the pRF
S = reshape(stimImage.im, nPixels, []);
fprintf('the dimensions of S are: %d (rows) by %d columns\n', size(S))


% the K matrix is the combination of stimulus and haemodynamics
K = H*S';
fprintf('the dimensions of H are: %d (rows) by %d columns\n', size(H))
fprintf('                  S'' are: %d (rows) by %d columns\n', size(S'))
fprintf('                  K = H*S'' are: %d (rows) by %d columns = H[%d by %d] * S''[%d by %d]\n', size(K),size(H),size(S'));

end
    
    