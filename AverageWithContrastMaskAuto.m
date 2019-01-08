%% Average angiogram volumes with applying contrast-based mask (auto)
% for each Y position, find the volume data (#=nvexc) where std/mean of z-MIP is smallest and reject it from averaging
% thr = relative number of B-scans that are artifacts
% namin = min of volume number to be averaged (when smaller than minv, average all volumes)
% nforcedavg = number of B-scans where avg number < namin so that all volumes were forced to be averaged

function DD = AverageWithContrastMaskAuto(DDD, nvexc, bFig)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jlee@optobrain.com');
end

if nargin < 4
    bFig = false;
end

    [nz,nx,ny,nv] = size(DDD);
    if nvexc >= nv
        error('nvexc is larger than nv.');
    end
    
    D = squeeze(max(DDD,[],1));  % MIP [nx ny nv]
    C = squeeze(std(D,1,1) ./ mean(D,1));  % [ny nv]
            if bFig
                figure;  plot(C);  xlabel('Y');  title('Contrast over different volumes');
                figure;  plot(min(C,[],2));  xlabel('Y');  title('Min Contrast');
            end
    
    DD = zeros(nz,nx,ny,'single');
    for iy=1:ny
        [~,is] = sort(C(iy,:));
        DD(:,:,iy) = mean(DDD(:,:,iy,is((nvexc+1):nv)),4);
    end
    
            