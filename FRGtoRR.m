%% Convert spectrum to reflectivity profile
% Required memory = FRG raw file X 4
% FRG0 [nk 1]

function RR = FRGtoRR(FRG, intpDk, bavgfrm, FRG0, window)

ConvertMMA2;

if nargin < 5
	window = ones(size(FRG,1),1);
end
if nargin < 4
	FRG0 = zeros(0);
end
if nargin < 3
	bavgfrm = false;
end
if nargin < 2
	intpDk = 0.25;
end

if size(FRG,3) > 1
    disp([datestr(now,'HH:MM') '  FRGtoRR started']);
    tic;
end

	%% Subtract mean

		[nk nx nf] = size(FRG);
		if numel(FRG0) > 0
			FRG = FRG ./ repmat(FRG0,[1 nx nf]) - 1;
		else
			if bavgfrm
				FRG = FRG ./ repmat(mean(FRG(:,:),2),[1 nx nf]) - 1;
			else
				for ifr=1:nf
					FRG(:,:,ifr) = FRG(:,:,ifr) ./ repmat(mean(FRG(:,:,ifr),2),[1 nx]) - 1;
				end
			end
        end

    %% multiply the window
        if nargin > 4
            FRG = FRG .* repmat(window,[1 nx nf]);
        end

	%% Interpolation of K
    %{

		if intpDk ~= 0
			k = linspace(1-intpDk/2, 1+intpDk/2, nk);
            lam = 1./k;
% 			lam = fliplr(lam);  % k should be decending :: when not, intpDk should be negative
			parfor ifr=1:nf
				FRG(:,:,ifr) = interp1(lam, FRG(:,:,ifr), linspace(min(lam),max(lam),length(lam)), 'spline');
% 				if (mod(ifr,ceil(nf/5)) == 0) && nf > 1
% 					disp([datestr(now,'HH:MM') '  FRGtoRR ' num2str(ifr) '/' num2str(nf) ' ...']);  
% 				end
			end			
		end
    %}


	%% Interpolation of K

		if intpDk ~= 0
            lam_range = intpDk;  % in percentage
            lam = linspace(1-lam_range/2,1+lam_range/2,nk);
            k = 1./lam;
            k_even = linspace(min(k),max(k),nk);
			parfor ifr=1:nf
				FRG(:,:,ifr) = interp1(k, FRG(:,:,ifr), k_even, 'spline');
% 				if (mod(ifr,ceil(nf/5)) == 0) && nf > 1
% 					disp([datestr(now,'HH:MM') '  FRGtoRR ' num2str(ifr) '/' num2str(nf) ' ...']);  
% 				end
			end			
		end


	%% Inverse Fourier transform
		
		nz = round(nk/2);
		RR = ones(nz,nx,nf,'single')*(1+1i);
		parfor ifr=1:nf
			RRy = ifft(FRG(:,:,ifr));
			RR(:,:,ifr) = RRy(1:nz,:);
        end
    

if size(FRG,3) > 1        
    disp([datestr(now,'HH:MM') '  FRGtoRR completed (' num2str(toc,2) ' s)']);
end


