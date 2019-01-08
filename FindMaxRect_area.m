function ret = FindMaxRect_area(B, c, nx, ny, nd)  % inverse of area for fminsearch

if nargin < 5
	nd = 0;
end

	c(1:2) = ceil(c(1:2));
	c(3:4) = floor(c(3:4));
	
	p1 = [0 0];  p2 = [0 0];
	p1x = min(max(c(1),1),nx/2-1);
	p1y = min(max(c(2),1),ny/2-1);
	p2x = min(max(c(3),nx/2+1),nx);
	p2y = min(max(c(4),ny/2+1),ny);
	
	b = B(p1x:p2x,p1y:p2y);
	
			if nd > 0
				subplot(4,4,nd+4);  cla;  colormap(gray);  hold on;
					PlotImage(B',0,[0 1],0);
					PlotBox(p1x:p2x,p1y:p2y,'r');
					title(mat2str([sum(b(:)) prod(B(:))]));
                pause(.1);
			end
	
	ret = 1/sum(b(:)) * prod(B(:));
	