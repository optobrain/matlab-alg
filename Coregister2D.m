%% find optimum rotation and translation to coregister a 2D map to a reference map

% fp = two featured point, fp(:,1) for translation and fp(:,2) for rotation

function [Ac, xcr, xcr0, r1, cor1, t1] = Coregister2D(A, A0, fp, fp0, id)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jlee@optobrain.com');
end

if nargin < 5
    id = 0;
end
	
	% initial guess
		t1 = fp(:,1) - fp0(:,1);  t1 = -t1;  % it works when this minus applied
		cor1 = fp(:,1);  % center of rotation
		r1 = 0;  % rotation angle
			v0 = fp0(:,2)-fp0(:,1);
			v1 = fp(:,2)-fp(:,1);
			r1 = acos(v1'*v0 / norm(v1) / norm(v0));  % positive when |r1| < pi/2
			if atan(v1) > atan(v0)  % it works when this opposite condition is tested
				r1 = -r1;
			end
			
	% optimize
		fitE2 = @(c) 1/Coregister2D_xcr(A, A0, c(1), cor1, [c(2) c(3)], id);
		fitC0 = [r1 t1(1) t1(2)];
		fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolX',0.1));
		r1 = fitC(1);  t1 = fitC(2:3);
	
	% coregistered A
		Ac = RotateImage(A,r1,cor1,t1);
		
	% xcr
		sqA = A(:);  sqA0 = A0(:);  sqAc = Ac(:);
		idx = find(~isnan(sqAc));
		xcr0 = GetXcr(sqA(idx),sqA0(idx));
		xcr = GetXcr(sqAc(idx),sqA0(idx));
		
				if id > 0
					B = ~isnan(Ac);
					cm = zeros(128,3);  cm(1:64,:) = gray(64);  cm(65:128,1) = linspace(0,1,64);
                    if id <= 2
                        jj = 2+id;
                    else
                        jj = 4+id;
                    end
					subplot(4,4,jj);  cla;  colormap(gca,cm);  hold on; 
						PlotImage(Rescale(A0',[0 0.9],[MeanLow(A0,0.1) MeanHigh(A0,0.9)]),0,[0 2],0);
						img = Rescale(Ac'.*B',[0.1 0.8],[MeanLow(sqAc(idx),0.1) MeanHigh(sqAc(idx),0.9)]);
						PlotImage(img+1,0,[0 2],0,[],img);
						set(gca,'CLim',[0 2]);
						title(['xcr : ' num2str(xcr0,2) ' => ' num2str(xcr,2)]);
				end
	
