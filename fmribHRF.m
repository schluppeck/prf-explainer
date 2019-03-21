% fmribHRF - hrf function used in by FMRIB / Oxford crew
%
%      usage: [ result ] = fmribHRF( t, p1, p2, p3, p4, p5, scaled )
%             [ result ] = fmribHRF( t, p, scaled )
%         by: denis schluppeck
%       date: 2006-05-01
%     inputs: t, p1, p2, p3, p4, p5 [or one vector for params], scaled [default=1]
%    outputs: p, result
%    purpose: calculate a haemodynamic response function according
%    to the equations given in the Jezzard et al book (page 252).
% 
%               H(t) = ( t / d1 ).^a1 .* exp( -(t-d1)./b1 ) 
%                      - c*( t / d2 ).^a2 .* exp( -(t-d2)./b2 ) 
% 
%   where di = ai*bi (default params a1, a2, b1, b2, c [6 12 0.9 0.9 0.35] )
%   
%       e.g.:  t = 0:30; % time in s
%              p = [6 12 0.9 0.9 0.35]
%              [result] = fmribHRF(t, p(1), p(2), p(3), p(4), p(5) )
%              [result] = fmribHRF(t, p )
% 
function [ result ] = fmribHRF(t, P, p2, p3, p4, p5, scaled)

if ~any( nargin == [6 7 2 3 1 0])
  help fmribHRF
  return
end

% all values defined.
if any(nargin == [6 7])
  p = zeros(5,1);
  p(1) = P; % that's it, all others are passed
  p(2) = p2; p(3) = p3; p(4) = p4; p(5) = p5;
  if nargin == 6
    scaled = 1;
  end % else scaled is given in input
end

% if a vector is passed.
if any(nargin == [2 3])
    p = P;
    if nargin == 2
      scaled = 1;
    else 
      % else scaled is given... but in p2
      scaled = p2;
    end
end


% if no parameters are passed use defaults from worsley chapter
% FMRIB book (p252)
if nargin < 2
    p = [6 12 0.9 0.9 0.35]; % a1, a2, b1, b2, c
    scaled = 1;
end

% if no time vector is given, specify 0:30 s
if nargin == 0
    t = 0:30;
end

% calculate the HRF
result = ( t ./ ( p(1)*p(3) )  ).^p(1) .* exp( -( t-p(1)*p(3) )./p(3) ) - ...
    p(5)*( t ./ ( p(2)*p(4) )  ).^p(2) .* exp( -( t-p(2)*p(4) )./p(4) );


if scaled
  % normalize to SUM=1
  result = result./sum(result);
end