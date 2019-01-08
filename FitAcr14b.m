%% ver. 14b :: force Vt value [scanning speed; m/s] with Vz=0
% bme: inital ms and mf from me
% bdv: D is determined by Ms and Mf during fitting  (3-coeff lsqcurvefit took really long time, 15x longer than ver. 14)
% V0: fixed V value [m/s]

function [Ms,Mf,D,V,A,R,GGf] = FitAcr14b(tau, GG, xzr, bme, bexp, Tol, vt0)

ConvertMMA2;

if (nargin < 7) vt0 = 0;  end;
if (nargin < 6) Tol = 1e-3;  end;
if (nargin < 5) bexp = true;  end;
if (nargin < 4) bme = true;  end;
if (nargin < 3) xzr = 1;  end;

	%% Prepare

		% constant
		k0 = 2*pi/1.31e-6;
%		dk = ( 2*pi/(1.31e-6-0.17e-6/2)-2*pi/(1.31e-6+0.17e-6/2) )/2/sqrt(2*log(2));
		dk = ( 2*pi/(1.31e-6-0.17e-6/2)-2*pi/(1.31e-6+0.17e-6/2) )/2*sqrt(2*log(2));
		n = 1.35;
		q = 2*n*k0;
		dq = 2*n*dk;
		h = dq/2;
		hx = h/xzr;
		dt = mean(diff(tau));

		% prepare
		[nz nx ny ntau] = size(GG);
		na = 3;
		Ms = zeros(nz,nx,ny,na);  Mf = Ms;  D = Ms;  V = Ms;  A = Ms;  R = zeros(nz,nx,ny,na);  GGf = ones(nz,nx,ny,ntau,na)*(1+1i);

	
	%% Fit

		t = reshape(tau(2:end),[ntau-1 1]);		
		tn = t / tau(end);						% tn(end) = 1 for TolX
		
		% guess the Ms0 as abs(g(end))
        Ms0 = abs(GG(:,:,:,end));

		% guess the initial Me0 with three points
		X = [tau(2)^2 tau(2) 1; tau(3)^2 tau(3) 1; tau(4)^2 tau(4) 1];
		C = Multiply(inv(X),shiftdim(abs(GG(:,:,:,2:4)),3),1);
		Me0 = 1 - min(max( shiftdim(C(3,:,:,:),1) ,0),1);
		clear X C;
		Mf0 = max(1-Ms0-Me0,0);

		% initial guess with FitAcr14bdv
		nN = 45;		% 0.1-0.9, 0.1-0.9
		Ms1 = zeros(1,nN);
		Mf1 = zeros(1,nN);
		iN = 0;
		for ms=1:9
			for mf=1:9
				if ms + mf <= 10
					iN = iN+1;
					Ms1(iN) = ms/10;
					Mf1(iN) = mf/10;
				end
			end
		end
		if ~bme
			[DD,RR] = FitAcr14bdv(tau, GG, Ms1, Mf1, q,h,hx,vt0);
			[r,IR] = max(RR,[],5);
		end

		for iz=1:nz
			for ix=1:nx
				for iy=1:ny

					% prepare
					g = squeeze(GG(iz,ix,iy,2:end));
					msi = zeros(3,1);  mfi = msi;  di = msi;  vti = msi;  vzi = msi;
					di0 = .5/q^2/tau(end);
					msf0 = 1-Me0(iz,ix,iy);

					% Guess 1: abs(g(end)) as initial value
					msi(1) = Ms0(iz,ix,iy);  mfi(1) = Mf0(iz,ix,iy);  di(1) = di0;

					% Guess 2: FitAcr14bdv as initial value
					if bme
						[DD,RR] = FitAcr14bdv(tau, GG(iz,ix,iy,:), Ms1*msf0, Mf1*msf0, q,h,hx,vt0);
						[r,ir] = max(RR,[],5);
						msi(2) = Ms1(ir)*msf0;  mfi(2) = Mf1(ir)*msf0;  di(2) = DD(1,1,1,1,ir);
					else
						msi(2) = Ms1(IR(iz,ix,iy));  mfi(2) = Mf1(IR(iz,ix,iy));  di(2) = DD(iz,ix,iy,1,IR(iz,ix,iy));
					end

					% general initial value
					msi(3) = msf0/2;  mfi(3) = msf0/2;  di(3) = di0;

					% fit
					for ja=1:3
                        Dmax = 100e-12;
                        if bexp
                            % this 2-coeff lsqcurvefit with external function took 36000 s (really; why?)
%                             fitF = @(c,tn) c(1) + min(c(2),1-c(1)) .* exp(-hx^2*vt0^2*t.^2-q^2.*squeeze(FitAcr14bdv(t,g,c(1),min(c(2),1-c(1)),q,h,hx,vt0)).*t);
                            %{
                            fitF = @(c,tn) c(1) + min(c(2),1-c(1)) .* exp(-hx^2*vt0^2*t.^2-q^2.*FitPolyn1( log(abs(g-c(1))./min(c(2),1-c(1)).*exp(hx^2*vt0^2*t.^2)) ,q^2*t,1,true).*t);
                            fitC0 = [msi(ja) mfi(ja)];
                            fitCL = [0 0];  fitCU = [1 1];
                            fitC = lsqcurvefit(fitF, fitC0, tn, g, fitCL, fitCU, optimset('Display','off','TolFun',Tol,'TolX',Tol));
                            ms = fitC(1);  mf = min(fitC(2),1-ms);  d = FitAcr14bdv(tau,GG(iz,ix,iy,:),ms,mf,q,h,hx,vt0);
                            %}
                            
                            fitE2 = @(c) sum( abs( min(abs(c(1)),1) + min(abs(c(2)),1-min(abs(c(1)),1)) .* exp(-hx^2*vt0^2*t.^2-q^2.*FitPolyn1( log(abs(g-c(1))./min(c(2),1-c(1)).*exp(hx^2*vt0^2*t.^2)) ,q^2*t,1,true).*t) - g ).^2 );
							fitC0 = [msi(ja) mfi(ja)];
							fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolFun',Tol,'TolX',Tol));
							ms = min(abs(fitC(1)),1);  mf = min(abs(fitC(2)),1-ms);  
                            d = FitAcr14bdv(tau,GG(iz,ix,iy,:),ms,mf,q,h,hx,vt0);
                            
                        else
                            % this 3-coeff lsqcurvefit took 900 s (ver. 14 took 50 s)
                            fitF = @(c,tn) c(1) + min(c(2),1-c(1)) .* exp(-hx^2*vt0^2*t.^2-c(3)*tn);
                            fitC0 = [msi(ja) mfi(ja) di(ja)];
                            fitCL = [0 0 0];  fitCU = [1 1 q^2*Dmax*tau(end)];
                            fitC = lsqcurvefit(fitF, fitC0, tn, g, fitCL, fitCU, optimset('Display','off','TolFun',Tol,'TolX',Tol));
                            ms = fitC(1);  mf = min(fitC(2),1-ms);  d = fitC(3)/q^2/tau(end);
                        end
                        v = vt0;
						
						Ms(iz,ix,iy,ja) = ms;  Mf(iz,ix,iy,ja) = mf;  D(iz,ix,iy,ja) = d;  V(iz,ix,iy,ja) = v;
						if (v > 0)  A(iz,ix,iy,ja) = vz/v;  end;
						GGf(iz,ix,iy,2:end,ja) = ms + mf*exp(-hx^2*vt0^2*t.^2-q^2*d*t);  GGf(iz,ix,iy,1,ja) = 1;
						R(iz,ix,iy,ja) = 1 - sum( abs(GG(iz,ix,iy,:)-GGf(iz,ix,iy,:,ja)).^2 ,4) / sum( abs(GG(iz,ix,iy,:)-mean(GG(iz,ix,iy,:),4)).^2 ,4);									
				
                    end
                    
				end
			end

			if (mod(iz,ceil(nz/10)) == 0)  
%				disp(['... FitAcr14 ' num2str(iz) '/' num2str(nz) '	' datestr(now,'HH:MM')]);  
			end
		end

%		GGf(:,:,:,1,:) = 1;
%		R = 1 - sum(abs(repmat(GG,[1 1 1 1 na])-GGf).^2,4) ./ repmat(sum(abs(GG-repmat(mean(GG,4),[1 1 1 size(GG,4)])).^2,4),[1 1 1 1 na]);
