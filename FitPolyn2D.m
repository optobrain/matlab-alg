%% Fit 2D data to f(x,y) = a*x + b*y + c or 2nd-order polynomial

% porder: polynomial order :: if 2 : x^2 + y^2 + x + y + c
% x and y: NOT a ROI, but the unit of x and y axis

function [f_fit, c] = FitPolyn2D(f, porder, x, y, mask, figno)

ConvertMMA2;

if nargin < 6
    figno = 0;
end
[nx,ny] = size(f);
if nargin < 5
    mask = ones(nx,ny);
end
if nargin < 4
    x = 1:nx;
    y = 1:ny;
end
if nargin < 2
	porder = 1;
end

    idx = find(mask(:));    
    [yy,xx] = meshgrid(y,x);
    xs = xx(:);  ys = yy(:);  fs = f(:);
	if porder == 3
        A = [xs(idx).^3 ys(idx).^3 xs(idx).^2 ys(idx).^2 xs(idx) ys(idx) ones(numel(idx),1)];
        c = inv(A'*A)*A'*fs(idx);
        f_fit = c(1)*xx.^3 + c(2)*yy.^3 + c(3)*xx.^2 + c(4)*yy.^2 + c(5)*xx + c(6)*yy + c(7);
    elseif porder == 2
        A = [xs(idx).^2 ys(idx).^2 xs(idx) ys(idx) ones(numel(idx),1)];
        c = inv(A'*A)*A'*fs(idx);
        f_fit = c(1)*xx.^2 + c(2)*yy.^2 + c(3)*xx + c(4)*yy + c(5);
    else
        A = [xs(idx) ys(idx) ones(numel(idx),1)];
        c = inv(A'*A)*A'*fs(idx);
        f_fit = c(1)*xx + c(2)*yy + c(3);
    end
    
    if figno > 0
        figure(figno);  clf;  colormap(jet);
            subplot(121);  hold on;  PlotImage(f',0,[.1 .9],1);  ax = gca;
            subplot(122);  hold on;  PlotImage(f_fit',0);  set(gca,'CLim',ax.CLim);
    end
