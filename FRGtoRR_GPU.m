%% Convert spectrum to reflectivity profile
% Required memory = FRG raw file X 4
% FRG0 [nk 1]

function RR = FRGtoRR_GPU(FRG, intpDk, bavgfrm, FRG0)

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

    gFRG = gpuArray(FRG);
    
	%% Subtract mean

		[nk nx nf] = size(gFRG);
		if numel(FRG0) > 0
			gFRG = gFRG - repmat(FRG0,[1 nx nf]);
		else
			if bavgfrm
				gFRG = gFRG - repmat(mean(gFRG(:,:),2),[1 nx nf]);
			else
				for ifr=1:nf
					gFRG(:,:,ifr) = gFRG(:,:,ifr) - repmat(mean(gFRG(:,:,ifr),2),[1 nx]);
				end
			end
		end


	%% Interpolation of K :: takes a long time

		if intpDk ~= 0
			k = linspace(1-intpDk/2, 1+intpDk/2, nk);
			lam = fliplr(1./k);  glam = gpuArray(lam);
            lamLin = linspace(lam(1),lam(end),length(lam));  glamLin = gpuArray(lamLin);
			for ifr=1:nf
				gFRG(:,:,ifr) = interp1(glam, gFRG(:,:,ifr), glamLin, 'linear');
			end			
		end


	%% Inverse Fourier transform
		
		nz = round(nk/2);
		RR = ones(nz,nx,nf,'single')*(1+1i);  RRy = ones(nk,nx,'single')*(1+1i);
        gRR = gpuArray(RR);  gRRy = gpuArray(RRy);
		for ifr=1:nf
			gRRy = ifft(gFRG(:,:,ifr));
			gRR(:,:,ifr) = gRRy(1:nz,:);
        end
        wait(gpu);
        RR = gather(gRR);
    

if size(FRG,3) > 1        
    disp([datestr(now,'HH:MM') '  FRGtoRR_GPU (' num2str(toc,2) ' s)']);
end


