%% find max-area rectangle within B=1 area (used to find common ROI from the coregistered MIPs), based on rotating rod

% oini = option of initial guess 1=half box, 2=diagonal intersection, 3=swipe angle [default=3]

% ret = two points of the rectangle. ret(:,1) ret(:,2)

function ret = FindMaxRect(B, oini, nd)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jonghwan_lee@brown.edu');
end

if nargin < 3
    nd = 0;
end
if nargin <2
	oini = 3;
end
	
	% initial guess
		[nx ny] = size(B);
		if oini == 1
			p1 = [nx ny]*1/4;
			p2 = [nx ny]*3/4;
			
		elseif oini == 2
			n = min(nx,ny);
			d = diag(B(1:n,1:n));
			idx = find(d);
			p1 = [1 1]*idx(1);
			d = diag(B(nx-n+1:nx,ny-n+1:ny));  d = d(end:-1:1);
			idx = find(d);
			p2 = [nx-idx(1)+1 ny-idx(1)+1];
			
		elseif oini == 3
			% find the angle of the rotating rod for max area
			[y,x] = meshgrid(1:ny,1:nx);
			cp = [mean(x(:).*B(:)) mean(y(:).*B(:))];  % center point
			
			x = 1:nx;
			ss = pi/2*[0.3:0.05:0.7];  ns = length(ss);  a = zeros(ns,1);
			for is=1:ns
				s = ss(is);
				y = tan(s)*(x-cp(1))+cp(2);  y = min(max(round(y),1),ny);
				b = zeros(nx,1);
				for ix=1:nx
					b(ix) = B(ix,y(ix));
				end
				idx = find(b);
				p1 = [idx(1) y(idx(1))];
				p2 = [idx(end) y(idx(end))];
				a(is) = prod(p2-p1);
			end
			[m im] = max(a);
			s = ss(im);

					if nd > 0
						subplot(4,4,nd+3);  cla;  plot(ss,a);  line(s,m,'marker','o');
					end

			% initial guess
			y = tan(s)*(x-cp(1))+cp(2);  y = min(max(round(y),1),ny);
			b = zeros(nx,1);
			for ix=1:nx
				b(ix) = B(ix,y(ix));
			end
			idx = find(b);
			p1 = [idx(1) y(idx(1))];
			p2 = [idx(end) y(idx(end))];
		end
			
	% optimize
		fitE2 = @(c) FindMaxRect_area(B,c,nx,ny,nd);
		fitC0 = [p1(1) p1(2) p2(1) p2(2)];
		fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolX',0.5));
		fitC(1:2) = ceil(fitC(1:2));  
		fitC(3:4) = floor(fitC(3:4));
		p1 = [ min(max(fitC(1)+1,1),nx/2) min(max(fitC(2)+2,1),ny/2) ];
		p2 = [ min(max(fitC(3)-2,nx/2),nx) min(max(fitC(4)-2,ny/2),ny) ];
		
	ret = zeros(2,2);
	ret(:,1) = p1';
	ret(:,2) = p2';
	
	
