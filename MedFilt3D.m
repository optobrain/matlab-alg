%% median filter 3D

% bpar = use parallel ?
% bwb = show waitbar?

function ret = MedFilt3D(I, kern, bpar)

ConvertMMA2;

if nargin < 3
	bpar = false;
end
if nargin < 2
	kern = [3 3 3];
end

	ret = I;
	edg = (kern-1)/2;
	[nz nx ny] = size(I);
%     disp('... MedFilt3D running');
    tic;
	if bpar
		parfor iz=1:nz
			for ix=1:nx
				for iy=1:ny
					img = I( min(max(iz+[-edg(1):edg(1)],1),nz) , min(max(ix+[-edg(2):edg(2)],1),nx) , min(max(iy+[-edg(3):edg(3)],1),ny) );
					ret(iz,ix,iy) = median(img(:));
				end
			end
		end
	else
		for iz=1:nz
			for ix=1:nx
				for iy=1:ny
					img = I( min(max(iz+[-edg(1):edg(1)],1),nz) , min(max(ix+[-edg(2):edg(2)],1),nx) , min(max(iy+[-edg(3):edg(3)],1),ny) );
					ret(iz,ix,iy) = median(img(:));
				end
			end
		end
    end
%     disp(['... MedFilt3D completed after ' num2str(round(toc/60)) ' min']);




