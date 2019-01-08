%% Fit 1-channel data with Gaussian function
% y = a exp( -(x-b)^2 / 2 / c^2 ) + d
% c0 = initial guess
% cmin, cmax
% bDzero = force the offset d=0?

function [a,b,c,d,R2,yfit] = FitGaussian(x, y, fitC0, fitCmin, fitCmax, bDzero)

ConvertMMA2;

if nargin < 6
    bDzero = false;
end

    if bDzero
        d = 0;
        fitF = @(c,x) c(1)*exp(-(x-c(2)).^2 /2/c(3)^2);
        if nargin < 5 || isempty(fitCmin)
            [fitC,resnorm] = lsqcurvefit(fitF,fitC0(1:3),x(:),y(:), optimset('Display','off'));
        else
            [fitC,resnorm] = lsqcurvefit(fitF,fitC0(1:3),x(:),y(:),fitCmin,fitCmax, optimset('Display','off'));
        end
        a = fitC(1);  b = fitC(2);  c = fitC(3);
    else
        fitF = @(c,x) c(1)*exp(-(x-c(2)).^2 /2/c(3)^2) + c(4);
        if nargin < 5 || isempty(fitCmin)
            [fitC,resnorm] = lsqcurvefit(fitF,fitC0,x(:),y(:), optimset('Display','off'));
        else
            [fitC,resnorm] = lsqcurvefit(fitF,fitC0,x(:),y(:),fitCmin,fitCmax, optimset('Display','off'));
        end
        a = fitC(1);  b = fitC(2);  c = fitC(3);  d = fitC(4);
    end
    yfit = fitF(fitC,x(:));
    R2 = 1 - sum( (y-yfit).^2 ) / sum( (y-mean(y)).^2 );
    
    