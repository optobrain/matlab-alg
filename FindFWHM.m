%% Find FWHM from given I(z) in terms of pixel
% om: option number: 1=simply by looking at Max/2 after upsampling by 16

function [fwhm,im1,im2] = FindFWHM(Iz, om)

ConvertMMA2;

if nargin < 2
	om = 1;
end

    if om == 1
        % find FWHM by looking for the half maximum
        zm = 1:length(Iz);
        zmu = [1:length(zm)*8]/8;
        Izu = interp1(zm,abs(Iz),zmu);
        [m im] = max(Izu);
        [m1 im1] = min(abs(Izu(1:im)-m/2));  im1 = im1/8;
        [m2 im2] = min(abs(Izu(im:end)-m/2));  im2 = (im2 + im-1)/8;
        fwhm = im2-im1;
    end
    
