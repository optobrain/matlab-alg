%% Calculate sharpness from the image
%% 03/30/10, Jonghwan, based on Vivek's code

function ret = GetSharpness(img, intpSharp)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jlee@optobrain.com');
end

if nargin < 2
	intpSharp = 2;		% default: Maximum entropy
end
	

	% Local variance
	if intpSharp == 1
		windowSize = 5;
		ret = 0;   
		for (igs=1:(size(img,1)-windowSize))
			switch (ndims(img))
				case (1)
					ret = ret + var( img(igs:(igs+windowSize)) );
				case (2)
					ret = ret + mean(var( img(igs:(igs+windowSize),:) ));
				case (3)
					ret = ret + mean(mean(var( img(igs:(igs+windowSize),:,:) )));
			end
		end


	% 1/|Entropy|
	elseif intpSharp == 2
		switch ndims(img)
			case (1)
				img = img.^2/sum(img.^2);
				ret = sum(img.*log(img+eps));
			case (2)
				img = img.^2/sum(sum(img.^2));
				ret = sum(sum(img.*log(img+eps)));  % entropy = < img * log(img) >_x,y
			case (3)
				img = img.^2/sum(sum(sum(img.^2)));
				ret = sum(sum(sum(img.*log(img+eps))));
		end
		ret = 1/abs(ret);
		
	% 1/Cross-correlation (how the cross-correlation as the function of lag is sharp)
	elseif intpSharp == 3
		maxlag = 11;
		corr = zeros(2*maxlag+1,1);
		switch (ndims(img))
			case (1)
				corr = xcorr(img,maxlag,'unbiased');
			case (2)
				for (ix=1:size(img,2))
					corr = corr + xcorr(img(:,ix),maxlag,'unbiased');
				end
			case (3)
				for (ix=1:size(img,2))
					for (iy=1:size(img,3))
						corr = corr + xcorr(img(:,ix,iy),maxlag,'unbiased');
					end
				end
		end
%		ret = sum(corr)/max(corr);
%		ret = 1/ret;
		ret = max(corr)/sum(corr);

	end

    