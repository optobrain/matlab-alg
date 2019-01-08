%% Get spectrum and f for given data and t

% Y : data
% t : time
% dim : dimension of time
% ftr : freq to time ratio = upsampling of frequency domain with zeropadding
% bnorm : normalize?
% bhamm : hamming window?


function [Yf f] = GetSpectrum(Y, t, dim, ftr, bnorm, bhamm)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jlee@optobrain.com');
end

if nargin < 6
	bhamm = 0;
end
if nargin < 5
	bnorm = 0;
end
if nargin < 4
	ftr = 1;
end
if nargin < 3
	dim = 1;
end

	% size
		nt = length(t);  ndim = ndims(Y);
		if size(Y,dim) ~= nt
			error('size(Y,dim) ~= length(t).');
		end
		sz = size(Y);								% [nz nx ny nt]
		sz1 = sz;  sz1(dim) = 1;					% [nz nx ny 1]
		szt = ones(size(sz));  szt(dim) = nt;		% [1 1 1 nt]
		szu = ones(size(sz));  szu(dim) = ftr;		% [1 1 1 ftr]
		sz2 = sz;  sz2(dim) = nt*ftr;				% [nz nx ny nt*ftr]

	% normalize
		if bnorm
			Y = ( Y - repmat(mean(Y,dim),szt) ) ./ repmat(std(Y,1,dim),szt);
		end

	% hamming
		if bhamm
			Y = Y .* repmat(reshape(hamming(nt),szt),sz1);
		end		

	% zeropadding
		if ftr > 1
			Y = repmat(Y,szu);  
			if dim > 1
				Y = shiftdim(Y,dim-1);
				Y(nt+1:end,:) = 0;
				Y = reshape(shiftdim(Y,ndim+1-dim),sz2);
			else
				Y(nt+1:end,:) = 0;
			end
		end
	
	% fft
		f = 1/mean(diff(t(:)))/2*linspace(0,1,nt*ftr/2+1);
		Yf = fft(Y,[],dim)/(nt/2);
