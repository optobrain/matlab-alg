%% Fix the focus curvature in an OCT volume data
% porder : polynomial order to fit the curvature (default: 2)
% binterp : do interpolation?
% Zf: Zf can be given from II and then used for fixing other volume data

function [IIfix, Zf] = FlattenFocusCurvature(II, Zf, porder, binterp, fno, xx, yy)

ConvertMMA2;

[nz,nx,ny] = size(II);
if nargin < 7
    xx = 1:nx;
    yy = 1:ny;
end
mask = zeros(nx,ny);
mask(xx,yy) = 1;
if nargin < 5
    fno = 0;
end
if nargin < 4
    binterp = true;
end
if nargin < 3
    porder = 2;
end
    
    iry = [2 round(ny/2) ny-1];  nyy = length(iry);
        if fno > 0
            figure(fno);  clf;  colormap(gray);
            for ii=1:nyy
                subplot(4,nyy,ii);  cla;  PlotImage(log(II(:,:,iry(ii))));  xlabel('x');  title(['iy=' num2str(iry(ii))]);
            end
            subplot(4,nyy,[1:nyy]+2*nyy);  cla;  plot(mean(II(:,:),2));
        end
            
    [M,Zm] = max(II,[],1);  M = squeeze(M);  Zm = squeeze(Zm);
        if fno > 0
            subplot(4,nyy,3*nyy+1);  cla;  colormap(gca,'jet');  hold on;  PlotImage(Zm');  PlotBox(xx,yy,'w',2);  xlabel('x');  ylabel('y');  ax = gca;  colorbar;
        end
    if nargin < 2 || numel(Zf) == 0
        Zf = FitPolyn2D(Zm,porder,1:nx,1:ny,mask); 
    end
        if fno > 0
            subplot(4,nyy,3*nyy+2);  cla;  colormap(gca,'jet');  hold on;  PlotImage(Zf');  title(['err=' num2str( mean(mean((Zm-Zf).^2))/mean(Zm(:)) )]);  set(gca,'CLim',ax.CLim);  colorbar;
        end
        
    Zfmin = min(Zf(:));
    if Zfmin < max(nz/10,20)  % the interpolation below selects from iz = Zf - Zfmin, so Zfmin is # of upper z's than the focus
        Zfmin = Zfmin + max(nz/10,20);
    end
    Zf = Zf - Zfmin;
    nzc = ceil(max(Zf(:)) - min(Zf(:)));
    IIfix = zeros(nz-nzc,nx,ny);
    parfor ix=1:nx
        for iy=1:ny
            if binterp
%                 IIfix(:,ix,iy) = interp1(1:nz,II(:,ix,iy),[1:nz-nzc]+Zf(ix,iy)-min(Zf(:)));
                IIfix(:,ix,iy) = interp1(1:nz,II(:,ix,iy),[1:nz-nzc]+Zf(ix,iy));
            else
%                 IIfix(:,ix,iy) = II([1:nz-nzc]+round(Zf(ix,iy)-min(Zf(:))),ix,iy);
                IIfix(:,ix,iy) = II([1:nz-nzc]+round(Zf(ix,iy)),ix,iy);
            end
        end
    end
        if fno > 0
            for ii=1:nyy
                subplot(4,nyy,ii+nyy);  PlotImage(log(IIfix(:,:,iry(ii))));  xlabel('x');
            end
            subplot(4,nyy,[1:nyy]+2*nyy);  line([1:size(IIfix,1)]+nzc,mean(IIfix(:,:),2),'color','r');  set(gca,'YScale','log');  set(gca,'YTick',10.^[1:10]);  axis tight;  grid on;
        end
        
    Zf = Zf + Zfmin;
