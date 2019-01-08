%% Average angiogram volumes with applying contrast-based mask
% find B-scans where std/mean of z-MIP is smaller than lower thr% and reject them from averaging
% thr = relative number of B-scans that are artifacts
% namin = min of volume number to be averaged (when smaller than minv, average all volumes)
% nforcedavg = number of B-scans where avg number < namin so that all volumes were forced to be averaged

function [AA, M, nforcedavg] = AverageWithContrastMask(AAA, thr, namin, bFig)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jlee@optobrain.com');
end

if nargin < 4
    bFig = false;
end

    [nz nx ny nv] = size(AAA);
    AAA = shiftdim(AAA,2);
    C = std(AAA(:,:,:),1,3) ./ mean(AAA(:,:,:),3);
    AAA = shiftdim(AAA,2);
            if bFig
                figure;  colormap(gray);  subplot(1,2,1);  hold on;  imagesc(C);  axis tight;  colorbar;  title('Contrast');  ylabel('Y');  xlabel('vol');
            end
    M = ones(ny,nv);
    for iv=1:nv
        C1 = sort(C(:,iv));
        M(C(:,iv)<C1(round(ny*thr)),iv) = 0;
    end
    nforcedavg = length(find(sum(M,2)<namin));
    M(sum(M,2)<namin,:) = 1;  % for iy for avgN < 4 :: better to average all volumes
    AA = zeros(nz,nx,ny);
    for iy=1:ny
        AA(:,:,iy) = sum( AAA(:,:,iy,:) .* repmat(reshape(M(iy,:),[1 1 1 nv]),[nz nx 1 1]) , 4) / sum(M(iy,:));
    end
            if bFig
                subplot(1,2,2);  hold on;  imagesc(M);  axis tight;  title(num2str(nforcedavg));
            end
            