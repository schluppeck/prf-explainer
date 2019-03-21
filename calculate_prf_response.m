function [out] = calculate_prf_response(stimImage,pRFhat, method, hrf)
%calculate_prf_response - dot product over many time frames (loop or reshaped)
%   
% ds 2019-03-14

% if no method is passed in, set default
if nargin < 3 || isempty(method)
    method = 'reshape';
end

if strcmpi(method, 'reshape')
    % idea: consider image as a 1d array!
    % first two dimenions: image size;, third: time
    nPixels = numel(stimImage.im(:,:,1));   
    
    % "vectorise" stimImage and also the pRF
    stimImage.imUnwrapped = reshape(stimImage.im, nPixels, []);
    pRFhatUnwrapped = pRFhat(:);
    
    % calculate response with one matrix product
    r = stimImage.imUnwrapped' * pRFhat(:);
  
elseif strcmpi(method, 'loop')
    % idea: consider image as a 1d array!
    nFrames = size(stimImage.im,3);
    
    % allocate space and loop over frames
    r = nan(nFrames, 1);
    t = zeros(size(stimImage.x));
    for iFrame = 1:nFrames
        t = stimImage.im(:,:,iFrame) .* pRFhat;
        r(iFrame) = sum(t(:));      
    end
end

% if the user asks for the response w/ hrf included, calculate it
% this will save time, if not!
if nargin < 4
    % no HRF specified... output r
    out = r;
else
    % convolve
    r_hrf = conv(r, hrf);
    % and keep the correct length of r_hrf
    r_hrf = r_hrf(1:numel(r));
    
    % and demean
    out = r_hrf - mean(r_hrf(:));
end


end

