%% Get Hessian matrix and their eigen values and eigen vectors, using the 1st/2nd order derivative of Gaussian kernel

% DD [nz nx ny] : angiogram usually: should be 2D at least
% wids : array of guess widths of tubular structures
% bconv : conv H with ones(3,3,3) before eigenvalue

% EE [nz nx ny nd nw] : nd eiven values per position in the order of the smallest absolute value, for each Gaussian width
% VV [nz nx ny nd nw] : direction normal vector (nd = ndim) = the eigen vector with the smallest absolute eigen value, for each Gaussian width
% HH [nz nx ny nd nd nw] : Hessian matrix, for each Gaussian width
	
function [EE,VV,HH] = GetHessian(DD, wids, bconv, bfig)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jonghwan_lee@brown.edu');
end

if nargin < 4
    bfig = false;
end
if nargin < 3
    bconv = false;
end
[nz,nx,ny] = size(DD);
nd = ndims(DD);


	% Gaussian kernels :: we shouldn't do normalization
		nw = length(wids);
		nkern = max(wids)*2;
		kern = zeros(2*nkern+1,2,nw);	% 1=first order, 2=second order derivative
		for iw=1:nw
			k = -nkern:nkern;  wid = wids(iw);
			kern(:,1,iw) = 1/sqrt(2*pi)/wid .* (-k/wid^2) .* exp(-k.^2/2/wid^2);
			kern(:,2,iw) = 1/sqrt(2*pi)/wid .* ( (-1/wid^2) + (-k/wid^2).^2 ) .* exp(-k.^2/2/wid^2);
		end			
			if bfig
				figure(1);  clf;
					for iw=1:nw
						k = -nkern:nkern;  wid = wids(iw);
						g = 1/sqrt(2*pi)/wid .* exp(-k.^2/2/wid^2);
						g1 = 1/sqrt(2*pi)/wid .* (-k/wid^2) .* exp(-k.^2/2/wid^2);
						g2 = 1/sqrt(2*pi)/wid .* ( (-1/wid^2) + (-k/wid^2).^2 ) .* exp(-k.^2/2/wid^2);
						subplot(2,2,1);  line(k,g*wid^2,'color','k');  title('G'); 	% later, will normalize eigen value for comparison across different sig's : see Eq. 21 in 1998 MIAN Sato
						subplot(2,2,2);  line(k,g1*wid^2,'color','k');  title('d/dx G');
						subplot(2,2,3);  line(k,g2*wid^2,'color','k');  title('d^2/dx^2 G');
					end
			end

	% Hessian matrix
		HH = zeros(nz,nx,ny,nd,nd,nw,'single');
		for iw=1:nw  % no parfor (it's using huge memory)
			for d1=1:nd
				for d2=1:nd
					if d1 == d2		% H = Gxx
						kern1 = kern(:,2,iw);
						if d1 == 1
							kern1 = reshape(kern1,[2*nkern+1 1 1]);
						elseif d1 == 2							
							kern1 = reshape(kern1,[1 2*nkern+1 1]);
						else
							kern1 = reshape(kern1,[1 1 2*nkern+1]);
						end
						HH(:,:,:,d1,d2,iw) = convn(DD,kern1,'same');
					else			% H = Gxz
						kern1 = kern(:,1,iw);
						if d1 == 1
							kern1 = reshape(kern1,[2*nkern+1 1 1]);
						elseif d1 == 2							
							kern1 = reshape(kern1,[1 2*nkern+1 1]);
						else
							kern1 = reshape(kern1,[1 1 2*nkern+1]);
						end
						HH(:,:,:,d1,d2,iw) = convn(DD,kern1,'same');
						if d2 == 1
							kern1 = reshape(kern1,[2*nkern+1 1 1]);
						elseif d2 == 2							
							kern1 = reshape(kern1,[1 2*nkern+1 1]);
						else
							kern1 = reshape(kern1,[1 1 2*nkern+1]);
						end
						HH(:,:,:,d1,d2,iw) = convn(HH(:,:,:,d1,d2,iw),kern1,'same');						
					end
				end
			end
			disp(['... Hessian ' num2str(iw) '/' num2str(nw) ' ' datestr(now,'HH:MM')]);
		end
		if bconv
			HH = convn(HH,ones(3,3,min(3,ny))/3/3/min(3,ny),'same');
        end

	% Eigen values and eigen vectors
    tic;
		EE = zeros(nz,nx,ny,nd,nw,'single');		% 3 eigen values for each position & wid
		VV = zeros(nz,nx,ny,nd,nw,'single');		% direction = eigen vector with the smallest absolute eigen values
		for iw=1:nw  % no parfor
			for iz=1:nz
				for ix=1:nx
					for iy=1:ny
                        h = squeeze(HH(iz,ix,iy,:,:,iw));
                        if sum(sum(isnan(h))) > 0 || sum(sum(isinf(h))) > 0
                            d = [1 0 0; 0 1 0; 0 0 1];
                            v = [1 0 0; 0 1 0; 0 0 1];
                        else
    						[v,d] = eig(h);
                        end
						d = diag(d);
						[s,is] = sort(abs(d));
						EE(iz,ix,iy,:,iw) = d(is);
						VV(iz,ix,iy,:,iw) = v(:,is(1));
					end
				end
			end
			disp(['... Eigenvalue ' num2str(iw) '/' num2str(nw) ' ' datestr(now,'HH:MM')]);
        end
    toc;

			if bfig
				iy = round(ny/2);
				figure(2);  clf;  colormap(jet);
					for iw=1:nw
						subplot(nw,10,(iw-1)*10+1);  hold on;  imagesc(DD);  axis image;  set(gca,'XTick',[]);  set(gca,'YTick',[]);
						for d1=1:2
							for d2=1:2
								subplot(nw,10,(iw-1)*10+(d1-1)*2+d2+1);  hold on;  imagesc(HH(:,:,iy,d1,d2,iw));  axis image;  set(gca,'CLim',[-1 1]*.25);  ylabel(mat2str([d1 d2]));  set(gca,'XTick',[]);  set(gca,'YTick',[]);
							end
						end
						subplot(nw,10,(iw-1)*10+6);  hold on;
							imagesc(EE(:,:,iy,1,iw));  axis image;  set(gca,'CLim',[-1 1]*2);  set(gca,'XTick',[]);  set(gca,'YTick',[]);  ylabel(num2str(MeanLow(abs(EE(:,:,:,1,iw)),0.1),2));
						subplot(nw,10,(iw-1)*10+7);  hold on;
							imagesc(EE(:,:,iy,2,iw));  axis image;  set(gca,'CLim',[-1 1]*2);  set(gca,'XTick',[]);  set(gca,'YTick',[]);  ylabel(num2str(MeanLow(EE(:,:,:,2,iw),0.1),2));
						for id=1:nd
							subplot(nw,10,(iw-1)*10+7+id);  hold on;
								imagesc(VV(:,:,iy,id,iw));  axis image;  set(gca,'CLim',[-1 1]);  set(gca,'XTick',[]);  set(gca,'YTick',[]);  ylabel(['V_' num2str(id)]);
						end
					end
			end

