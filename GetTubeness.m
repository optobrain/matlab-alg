%% Get Tubeness with widths using the Hessian matrix

% DD [nz nx ny] : angiogram usually: should be 2D at least
% wids : array of guess widths of tubular structures
% para [alp gam12 gam23] : parameters for tubeness
% bconv : conv H & D with ones(3,3,3)

% TT [nz nx ny] : Tubeness = maximum D123 over wids in 1998 MIAN Sato)
% WW [nz nx ny] : The closest Gaussian width with the maximum T
% VV [nz nx ny nd] : direction normal vector (nd = ndim) = the eigen vector with the smallest absolute eigen value and the closest Gaussian width
	
function [TT,WW,VV] = GetTubeness(DD, wids, para, bconv, bfig)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jonghwan_lee@brown.edu');
end

if nargin < 5
    bfig = false;
end
if nargin < 4
    bconv = false;
end
if nargin < 3
    para = [.25 .5 .5];
end

% 		[nz,nx,ny] = size(DD);
% 		nd = ndims(DD);
		alp = para(1);  gam12 = para(2);  gam23 = para(3);

	% Hessian matrix and eigenvalues
		[EE,VV] = GetHessian(DD, wids, bconv, bfig);
		[nz,nx,ny,nd,nw] = size(EE);
		if bconv
			EE = convn(EE,ones(3,3,min(3,ny))/3/3/min(3,ny),'same');
		end
		Dn = EE .* repmat(reshape(wids.^2,[1 1 1 1 nw]),[nz nx ny nd 1]);	% normalize eigen value for comparison across different sig's : see Eq. 21 in 1998 MIAN Sato

	% D123 : Eq. 21 in 1998 MIAN Sato
		D123 = zeros(nz,nx,ny,nw,'single');  D123 = D123(:);
		sqD1 = Dn(:,:,:,1,:);  sqD1 = sqD1(:);
		sqD2 = Dn(:,:,:,2,:);  sqD2 = sqD2(:);
		if nd == 3
			sqD3 = Dn(:,:,:,3,:);  sqD3 = sqD3(:);

			idx = find(sqD3 < sqD2 & sqD2 < sqD1 & sqD1 <= 0);
			D123(idx) = abs(sqD3(idx)) .* (sqD2(idx)./sqD3(idx)).^gam23 .* (1+sqD1(idx)./abs(sqD2(idx))).^gam12;

			idx = find(sqD3 < sqD2 & sqD2 < 0 & 0 < sqD1 & sqD1 < abs(sqD2)/alp);
			D123(idx) = abs(sqD3(idx)) .* (sqD2(idx)./sqD3(idx)).^gam23 .* (1-alp*sqD1(idx)./abs(sqD2(idx))).^gam12;

		else	% lam3 = lam2 for 2D

			idx = find(sqD2 < sqD1 & sqD1 <= 0);
			D123(idx) = abs(sqD2(idx)) .* (1+sqD1(idx)./abs(sqD2(idx))).^gam12;

			idx = find(sqD2 < 0 & 0 < sqD1 & sqD1 < abs(sqD2)/alp);
			D123(idx) = abs(sqD2(idx)) .* (1-alp*sqD1(idx)./abs(sqD2(idx))).^gam12;
		end

		D123 = reshape(D123,[nz nx ny nw]);
				if bfig
					figure(1);  clf;  colormap(gray);
						for iw=1:nw
							subplot(1,nw,iw);  hold on;  img = max(D123(:,:,:,iw),[],3);  imagesc(img);  axis image;  set(gca,'XTick',[]);  set(gca,'YTick',[]);  colorbar;
						end
				end
	
	% T and W
		[TT,im] = max(D123,[],4);  
		WW = single(wids(im));
		Vm = zeros(nz,nx,ny,nd,'single');
		for iz=1:nz
			for ix=1:nx
				for iy=1:ny
					Vm(iz,ix,iy,:) = VV(iz,ix,iy,:,im(iz,ix,iy));
				end
			end
		end
		VV = Vm;


