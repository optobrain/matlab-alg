
function ret = CorrPhase(RR, it0, ncorrect, z, K, bFig)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jonghwan_lee@brown.edu');
end

if nargin < 6
	bFig = false;
end
if nargin < 5
	K = ones(size(RR,1),size(RR,2));
end
if nargin < 4
	z = 1:size(RR,1);
end
if nargin < 3
	ncorrect = 3;
end
if nargin < 2
	it0 = 1;
end

		% K should be > 0
		K = max(K,eps);

		% phase(zxzxz)
		for nphasecorr = 1:ncorrect
			% phase(z)
			I1 = mean( RR(z,:,:,:) .* repmat(conj(RR(z,:,:,it0)) .* K(z,:).^2 ,[1 1 1 size(RR,4)]) ,1);
			I1 = I1 ./ abs(I1);
			RR = RR ./ repmat(I1,[size(RR,1) 1 1 1]);
				% Plot
				if bFig
					I11 = I1(1,:,:,:);
					I1f = abs(fft(I11,[],4));
					figure;
						subplot(3,1,1);		imagesc(squeeze( real(I11) ));		colorbar;		title('real');
						subplot(3,1,2);		imagesc(squeeze( imag(I11) ));		colorbar;		title('imag');
						subplot(3,1,3);		imagesc(squeeze( log10(I1f) ));		colorbar;;		title('fft');
				end

%			if ncorrect > 1 && nphasecorr == ncorrect
			if nphasecorr == ncorrect
				break;
			end

			% phase(x)
			I1 = mean( RR .* repmat(conj(RR(:,:,:,it0)) .* K.^2 ,[1 1 1 size(RR,4)]) ,2);
			I1 = I1 ./ abs(I1);
			RR = RR ./ repmat(I1,[1 size(RR,2) 1 1]);
				% Plot
				if bFig
					I11 = I1(1,:,:,:);
					I1f = abs(fft(I11,[],4));
					figure;
						subplot(3,1,1);		imagesc(squeeze( real(I11) ));		colorbar;		title('real');
						subplot(3,1,2);		imagesc(squeeze( imag(I11) ));		colorbar;		title('imag');
						subplot(3,1,3);		imagesc(squeeze( log10(I1f) ));		colorbar;;		title('fft');
				end

			if ncorrect == 1
				break;
			end

		end

	ret = RR;

		
