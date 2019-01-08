%% Find COR Ro s.t. std of distances from Ro is minimized.
% If we minimize std/mean, then it would result in a very far point, because the farer point will lead to a smaller std/mean.

function ret = FindCOR(RR, bFig)

ConvertMMA2;

if nargin < 2
	bFig = false;
end

% tic;

	a = real(RR);
	b = imag(RR);

	A1 = std(a,1,4).^2;
	A2 = mean( (a - repmat(mean(a,4),[1 1 1 size(a,4) 1])) .* b ,4);
	A3 = mean( (b - repmat(mean(b,4),[1 1 1 size(b,4) 1])) .* a ,4);
	A4 = std(b,1,4).^2;
	
	B1 = 1/2 * ( mean(a.^3,4) - mean(a,4).*mean(a.^2,4) + mean( (a - repmat(mean(a,4),[1 1 1 size(a,4) 1])) .* b.^2 ,4) );
	B2 = 1/2 * ( mean(b.^3,4) - mean(b,4).*mean(b.^2,4) + mean( (b - repmat(mean(b,4),[1 1 1 size(b,4) 1])) .* a.^2 ,4) );

	COR = ( (A4.*B1 - A2.*B2) + 1i*(A1.*B2 - A3.*B1) ) ./ (A1.*A4 - A2.*A3);
    
    % COR is NaN when A1.*A4 - A2.*A3 = 0
    sz = size(COR);
    COR = COR(:);
    COR(isnan(COR)) = 0;
    ret = reshape(COR,sz);


		if bFig
				figure;
					Rf0 = mean(abs(RR(:,:,:,:,1) - repmat(ret,[1 1 1 size(RR,4) 1])),4);
					Rfstd = std(abs(RR(:,:,:,:,1) - repmat(ret,[1 1 1 size(RR,4) 1])),1,4);
					for im=1:4
						subplot(2,2,im);
							switch im
								case 1
									[m izx] = FindMax(Rf0);
								case 2
									[m izx] = FindMax(Rf0,'min');
								case 3
									[m izx] = FindMax(Rfstd);
								case 4
									[m izx] = FindMax(Rfstd,'min');
							end
							plot(ret(izx(1),izx(2),:,:,1),'Marker','o','color','k');
							PlotTrace( squeeze( Rf0(izx(1),izx(2),:,:,1) * exp(1i*linspace(1,4*pi,64)) + ret(izx(1),izx(2),:,:,1) ),'gray','.');
							PlotTrace( squeeze( RR(izx(1),izx(2),:,:,1) ),'jet','.');
							axis image;		set(gca,'XTick',[-1:1]);	set(gca,'YTick',[-1:1]);		grid on;
					end

        end

% disp([datestr(now,'HH:MM') '  FindCOR completed (' num2str(toc,2) ' s)']);

    