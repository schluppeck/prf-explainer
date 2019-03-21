function m = mgauss(p,x,y)
% mgauss - 2d gaussian (x,y, with params p x0,y0, sdxy)
%
%    inputs: p (parameter vector, [x0,y0, sdxy] )
%            x,y (meshgrid format)
%            
%           NB! params first, to make optimisation easier
%
% ds 2019-03-20 (break out of mrTools/mgl version)

% unpack parameters for readability
x0 = p(1);
y0 = p(2);
sdxy = p(3);

% compute 2d gaussian
m = exp(-(((x-x0).^2)/(2*(sdxy^2))+((y-y0).^2)/(2*(sdxy^2))));

% clamp small values to 0
m(m(:)<0.01) = 0;

end