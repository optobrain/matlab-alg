%% Rotate 2D image then translation

% I [nx ny]
% r : rotating angle (CCW)
% cor : center of rotation [ix iy]
% t : translation [tx ty]

function ret = RotateImage(I, r, cor, t)

ConvertMMA2;

if nargin < 4
    t = [0 0];
end
[nx ny] = size(I);
if nargin < 3 || numel(cor) == 0
    cor = [nx/2 ny/2];
end
    
    t = t(end:-1:1);  cor = cor(end:-1:1);  % i don't know why but this works correctly

    [x,y] = meshgrid(1:ny,1:nx);  % change this to [y,x] = meshgrid(1:ny,1:nx) then the above may not need.
    xi = x - cor(1);  yi = y - cor(2); 
    xi1 = xi*cos(r) - yi*sin(r);
    yi1 = xi*sin(r) + yi*cos(r);
    xi = xi1 + cor(1);  yi = yi1 + cor(2);

    xi = xi - t(1);  yi = yi - t(2);
    
    ret = interp2(x,y,I,xi,yi);  % don't use 'spline' because it did not return NaN
