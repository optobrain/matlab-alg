%% Find Z (depth map of MIP) with intensity mask

% porder: polynomial order :: if 2 : x^2 + y^2 + x + y + c
% ralv: relative area of large (surface) vessels for the mask

function ret = GetSurfaceTilt(D, porder, ralv, bFig)

ConvertMMA2;

if nargin < 4
    bFig = 0;
end
if nargin < 3
    ralv = 0.15;
end
if nargin < 2
	porder = 1;
end

	[nz nx ny] = size(D);
	[m z] = max(D,[],1);
%            figure(1);  clf;  colormap(jet);  hold on;  PlotImage(squeeze(m)',0);  colorbar;
%            figure(2);  clf;  colormap(jet);  hold on;  PlotImage(squeeze(z)',0);  colorbar;
    thr = sort(m(:),'descend');  thr = thr(round(numel(thr)*ralv));  mm = m>thr;
%            figure(3);  clf;  colormap(gray);  hold on;  PlotImage(squeeze(mm)',0);
%            figure(3);  clf;  colormap(jet);  hold on;  PlotImage(squeeze(z)',0,[0 200],0,[],squeeze(mm)');
    idx = find(mm(:));    
	if porder == 2
		if ny > 1
			[y x] = meshgrid(1:ny,1:nx);
			A = [x(idx).^2 y(idx).^2 x(idx) y(idx) ones(numel(idx),1)];
			c = inv(A'*A)*A'*z(idx);
			zm = c(1)*x.^2 + c(2)*y.^2 + c(3)*x + c(4)*y + c(5);
		else
			x = 1:nx;
			A = [x(idx).^2 x(idx) ones(numel(idx),1)];
			c = inv(A'*A)*A'*z(idx);
			zm = c(1)*x'.^2 + c(2)*x + c(3);
		end
	else
		if ny > 1
			[y x] = meshgrid(1:ny,1:nx);
			x = x(:);  y = y(:);  z = z(:);
			x = x(idx);  y = y(idx);  z = z(idx);
			c = inv([sum(x.^2) sum(x.*y) sum(x); sum(x.*y) sum(y.^2) sum(y); sum(x) sum(y) length(x)]) * [sum(z.*x) sum(z.*y) sum(z)]';
			[y x] = meshgrid(1:ny,1:nx);
			zm = c(1)*x + c(2)*y + c(3);
		else
			x = 1:nx;
            x = x(idx);
			c = inv([sum(x.^2) sum(x); sum(x) nx]) * [sum(z.*x) sum(z)]';
			x = 1:nx;
			zm = c(1)*x' + c(2);
		end
    end
    
    if bFig
        figure(11);  clf;  colormap(jet);
            subplot(221);  hold on;  PlotImage(squeeze(log(m))',0);
            subplot(222);  hold on;  PlotImage(squeeze(z)',0,[0 nz/2]);
            subplot(223);  hold on;  PlotImage(squeeze(z)',0,[0 nz/2],0,[],squeeze(mm)');
            subplot(224);  hold on;  PlotImage(squeeze(zm)',0,[0 nz/2]);
    end

ret = zm;
