%% extract four profiles from 2D img and fit the mean of four into normalized Gaussian
% d2 = (distance from the center)^2
% (1-b)*exp(-a*d2) + b

function [gb,gb0,d2,a,b,sig] = FitCrossSection(img)

ConvertMMA2;

    ned = floor((min(size(img))-1)/2);

    gb = zeros(ned,4);  gb(:,1) = img(ned+1,ned+1:end-1);  gb(:,2) = img(ned+1,ned+1:-1:2);  gb(:,3) = img(ned+1:-1:2,ned+1);  gb(:,4) = img(ned+1:end-1,ned+1);  
    gb = mean(gb,2);  gb0 = gb(1);  gb = gb/gb0;
    d2 = [0:ned-1]'.^2;  
        fitE2 = @(c) sum( ( (1-abs(c(2))).*exp(-abs(c(1)).*d2) + abs(c(2)) - gb ).^2 );
        fitC0 = [1/2/3^2 gb(end)];
        fitC = fminsearch(fitE2, fitC0, optimset('Display','off'));
        a = abs(fitC(1));  b = abs(fitC(2));
    sig = 1/sqrt(2*a);
