%% Fit 1-channel data with 2D Gaussian function
% z = a exp( -((x-b)^2 + (y-c)^2) / 2 / d^2 ) + e  [nx ny]
% [gy,gx] = meshgrid(1:ny,1:nx)
% c0 = initial guess
% bEzero = force the offset e=0?

function [a,b,c,d,e,R2,zfit] = FitGaussian2D(gx, gy, zdata, fitC0, bEzero)

ConvertMMA2;

if nargin < 5
    bEzero = false;
end

    if bEzero
        fitF = @(c,x,y) c(1)*exp(-((x-c(2)).^2+(y-c(3)).^2) /2/c(4)^2);
        fitE2 = @(c) mean(mean( ( fitF(c,gx,gy) - zdata ).^2 ));
        fitC = fminsearch(fitE2, fitC0, optimset('Display','off'));
        a = fitC(1);  b = fitC(2);  c = fitC(3);  d = fitC(4);  e = 0;
    else
        fitF = @(c,x,y) c(1)*exp(-((x-c(2)).^2+(y-c(3)).^2) /2/c(4)^2) + c(5);
        fitE2 = @(c) mean(mean( ( fitF(c,gx,gy) - zdata ).^2 ));
%         fitE2 = @(c) mean(mean( ( ( c(1)*exp(-((gx-c(2)).^2 + (gy-c(3)).^2) /2/c(4)^2) + c(5) ) - zdata ).^2 ));
        fitC = fminsearch(fitE2, fitC0, optimset('Display','off'));
        a = fitC(1);  b = fitC(2);  c = fitC(3);  d = fitC(4);  e = fitC(5);  
    end
    zfit = fitF([a b c d e],gx,gy);
    R2 = 1 - sum(sum( (zdata-zfit).^2 )) / sum(sum( (zdata-mean(zdata(:))).^2 ));
    
    