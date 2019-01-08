%% Convert spectrum to reflectivity profile
% Required memory = FRG raw file X 4
% FRG0 [nk 1]

function RR = FRGtoRR_GPU2(FRG, intpDk, bavgfrm, FRG0)

ConvertMMA2;

if nargin < 4
	FRG0 = zeros(0);
end
if nargin < 3
	bavgfrm = false;
end
if nargin < 2
	intpDk = 0.37;
end

if size(FRG,3) > 1
    disp([datestr(now,'HH:MM') '  FRGtoRR_GPU started']);
    tic;
end

    gpu = gpuDevice();

	%% Subtract mean

		[nk nx nf] = size(FRG);
		if numel(FRG0) > 0
			FRG = FRG - repmat(FRG0,[1 nx nf]);
		else
			if bavgfrm
				FRG = FRG - repmat(mean(FRG(:,:),2),[1 nx nf]);
			else
				for ifr=1:nf
					FRG(:,:,ifr) = FRG(:,:,ifr) - repmat(mean(FRG(:,:,ifr),2),[1 nx]);
				end
			end
		end


	%% Interpolation of K :: takes a long time

		if intpDk ~= 0
			k = linspace(1-intpDk/2, 1+intpDk/2, nk);
			lam = fliplr(1./k);  glam = gpuArray(lam);
            lamLin = linspace(lam(1),lam(end),length(lam));  glamLin = gpuArray(lamLin);
			for ifr=1:nf
                gFRG = gpuArray(FRG(:,:,ifr));
                gFRG = interp1(glam, gFRG, glamLin, 'linear');  % only linear is possible on GPU
%                 wait(gpu);
                FRG(:,:,ifr) = gather(gFRG);
			end			
		end


	%% Inverse Fourier transform
%{		
		nz = round(nk/2);
		RR = ones(nz,nx,nf,'single')*(1+1i);  RRy = ones(nk,nx,'single')*(1+1i);
        gRR = gpuArray(RR);  gRRy = gpuArray(RRy);
		for ifr=1:nf
			gRRy = ifft(gFRG(:,:,ifr));
			gRR(:,:,ifr) = gRRy(1:nz,:);
        end
        wait(gpu);
        RR = gather(gRR);
%}    
		nz = round(nk/2);
		RR = ones(nz,nx,nf,'single')*(1+1i);
		parfor ifr=1:nf
			RRy = ifft(FRG(:,:,ifr));
			RR(:,:,ifr) = RRy(1:nz,:);
        end

if size(FRG,3) > 1        
    disp([datestr(now,'HH:MM') '  FRGtoRR_GPU (' num2str(toc,2) ' s)']);
end


