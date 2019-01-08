% bme: inital ms and mf from me
% bexp: fit ms and mf only using FitAcr14exp
% bw: weighting radius and exponential?

function [Ms Mf D V A R GGf] = FitAcr15(tau, GG, xzr, bme, bexp, bw, Tol)

ConvertMMA2;

if (nargin < 7) Tol = 1e-3;  end;
if (nargin < 6) bw = 1;  end;
if (nargin < 5) bexp = 1;  end;
if (nargin < 4) bme = 1;  end;
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
		
		% guess the Ms0 (COR) and Vz0
		Ms0 = min(max( real( FindCOR(GG(:,:,:,ceil(end/2):end)) ) ,0),1);		% real was better than abs
		GGa = unwrap(angle(GG-repmat(Ms0,[1 1 1 ntau])),[],4);			
		T = repmat(reshape(tau,[1 1 1 ntau]),[nz nx ny 1]);
		if bw == 1
			Vz0 = sum(GGa.*T.*abs(GG-repmat(Ms0,[1 1 1 ntau])),4) ./ sum(T.^2.*abs(GG-repmat(Ms0,[1 1 1 ntau])),4) /q;
		else
			Vz0 = sum(GGa.*T,4) ./ sum(T.^2,4) /q;		% = Multiply(GGa,tau') / (tau*tau') /q
		end
		clear GGa T;

		% guess the initial Me0 with three points
		X = [tau(2)^2 tau(2) 1; tau(3)^2 tau(3) 1; tau(4)^2 tau(4) 1];
		C = Multiply(inv(X),shiftdim(abs(GG(:,:,:,2:4)),3),1);
		Me0 = 1 - min(max( shiftdim(C(3,:,:,:),1) ,0),1);
		clear X C;
		Mf0 = max(1-Ms0-Me0,0);

		% initial guess with FitAcr14dv
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
		if bme == 0
			[DD VT VZ RR] = FitAcr14dv(tau, GG, Ms1, Mf1, q,h,hx);
			[r IR] = max(RR,[],5);
		end

		for iz=1:nz
			for ix=1:nx
				for iy=1:ny

					% prepare
					g = squeeze(GG(iz,ix,iy,2:end));
					msi = zeros(3,1);  mfi = msi;  di = msi;  vti = msi;  vzi = msi;
					di0 = .5/q^2/tau(end);
					vi0 = .5/h/tau(end);
					msf0 = 1-Me0(iz,ix,iy);

					% COR as initial value
					msi(1) = Ms0(iz,ix,iy);  mfi(1) = Mf0(iz,ix,iy);  di(1) = di0;  vti(1) = Vz0(iz,ix,iy);  vzi(1) = Vz0(iz,ix,iy);

					% FitAcr14dv as initial value
					if bme == 1
						[DD VT VZ RR] = FitAcr14dv(tau, GG(iz,ix,iy,:), Ms1*msf0, Mf1*msf0, q,h,hx);
						[r ir] = max(RR,[],5);
						msi(2) = Ms1(ir)*msf0;  mfi(2) = Mf1(ir)*msf0;  di(2) = DD(1,1,1,1,ir);  vti(2) = VT(1,1,1,1,ir);  vzi(2) = VZ(1,1,1,ir);
					else
						msi(2) = Ms1(IR(iz,ix,iy));  mfi(2) = Mf1(IR(iz,ix,iy));  di(2) = DD(iz,ix,iy,1,IR(iz,ix,iy));  vti(2) = VT(iz,ix,iy,1,IR(iz,ix,iy));  vzi(2) = VZ(iz,ix,iy,IR(iz,ix,iy));
					end

					% general initial value
					msi(3) = msf0/2;  mfi(3) = msf0/2;  di(3) = di0;  vti(3) = vi0/sqrt(2);  vzi(3) = vi0/sqrt(2);

					% fit
					for ja=1:3
						if bexp == 1		% use FitAcr14exp
							fitE2 = @(c) sum( abs( min(abs(c(1)),1) + min(abs(c(2)),1-min(abs(c(1)),1)) .* FitAcr14exp(tn, g, min(abs(c(1)),1), min(abs(c(2)),1-min(abs(c(1)),1)), q,h,bw) - g ).^2 );
% 							fitE2 = @(c) sum( abs( min(abs(c(1)),1-min(abs(c(2)),1)) + min(abs(c(2)),1) .* FitAcr14exp(tn, g, min(abs(c(1)),1-min(abs(c(2)),1)), min(abs(c(2)),1), q,h,bw) - g ).^2 );  % same result
							fitC0 = [msi(ja) mfi(ja)];
							fitE2 = @(c) sum( abs( min(abs(0.527),1) + min(abs(c(1)),1-min(abs(0.527),1)) .* FitAcr14exp(tn, g, min(abs(0.527),1), min(abs(c(1)),1-min(abs(0.527),1)), q,h,bw) - g ).^2 );
% 							fitE2 = @(c) sum( abs( min(abs(c(1)),1-min(abs(c(2)),1)) + min(abs(c(2)),1) .* FitAcr14exp(tn, g, min(abs(c(1)),1-min(abs(c(2)),1)), min(abs(c(2)),1), q,h,bw) - g ).^2 );  % same result
							fitC0 = [mfi(ja)];
							fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolFun',Tol,'TolX',Tol));
% 							fun = @(c,tn) c(1) + min(c(2),1-c(1)) .* FitAcr14exp(tn, g, c(1), min(c(2),1-c(1)), q,h,bw);
% 							fitC0 = [msi(ja) mfi(ja)];
% 							fitC = lsqcurvefit(fun, fitC0, tn, g, [0 0], [1 1], optimset('Display','off','TolFun',Tol,'TolX',Tol));
%                             using lsqcurvefit resulted in lower R2
							ms = min(abs(0.527),1);  mf = min(abs(fitC(1)),1-ms);  
							[Exp d vt vz] = FitAcr14exp(tn, g, ms, mf, q,h,bw);
							d = d/q^2/tau(end);  vt = vt/hx/tau(end);  vz = vz/q/tau(end);  v = sqrt(vt^2+vz^2);
						else
							if xzr == 1
								fitE2 = @(c) sum( abs( min(abs(c(1)),1) + min(abs(c(2)),1-min(abs(c(1)),1)).*exp( -c(3).^2.*tn.^2 -abs(c(4)).*tn ).*exp(1i*q* sum(unwrap(angle(g-min(abs(c(1)),1))).*t)/sum(t.^2)/q *t) - g ).^2 );
								fitC0 = [msi(ja) mfi(ja) h*sqrt(vti(ja)^2+vzi(ja)^2)*tau(end) q^2*di(ja)*tau(end)];
							else
								fitE2 = @(c) sum( abs( min(abs(c(1)),1) + min(abs(c(2)),1-min(abs(c(1)),1)).*exp( -h^2*(sum(unwrap(angle(g-min(abs(c(1)),1))).*t)/sum(t.^2)/q)^2*t.^2 -c(3).^2.*tn.^2 -abs(c(4)).*tn ).*exp(1i*q* sum(unwrap(angle(g-min(abs(c(1)),1))).*t)/sum(t.^2)/q *t) - g ).^2 );
								fitC0 = [msi(ja) mfi(ja) hx*vti(ja)*tau(end) q^2*di(ja)*tau(end)];
							end
							fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolFun',Tol,'TolX',Tol));
% 							if xzr == 1
% 								fun = @(c,tn) c(1) + min(c(2),1-c(1)).*exp( -c(3).^2.*tn.^2 -c(4).*tn ).*exp(1i*q* sum(unwrap(angle(g-c(1))).*t)/sum(t.^2)/q *t);
% 								fitC0 = [msi(ja) mfi(ja) h*sqrt(vti(ja)^2+vzi(ja)^2)*tau(end) q^2*di(ja)*tau(end)];
% 							else
% 								fitE2 = @(c) sum( abs( min(abs(c(1)),1) + min(abs(c(2)),1-min(abs(c(1)),1)).*exp( -h^2*(sum(unwrap(angle(g-min(abs(c(1)),1))).*t)/sum(t.^2)/q)^2*t.^2 -c(3).^2.*tn.^2 -abs(c(4)).*tn ).*exp(1i*q* sum(unwrap(angle(g-min(abs(c(1)),1))).*t)/sum(t.^2)/q *t) - g ).^2 );
% 								fitC0 = [msi(ja) mfi(ja) hx*vti(ja)*tau(end) q^2*di(ja)*tau(end)];
% 							end
% 							fitC = lsqcurvefit(fun, fitC0, tn, g, [0 0 0 0], [1 1 h*sqrt((10e-3)^2+(10e-3)^2)*tau(end) q^2*20e-12*tau(end)], optimset('Display','off','TolFun',Tol,'TolX',Tol));
							ms = min(abs(fitC(1)),1);  mf = min(abs(fitC(2)),1-ms);  d = abs(fitC(4))/q^2/tau(end);
							vz = sum(unwrap(angle(g-ms)).*t)/sum(t.^2)/q;
							if xzr == 1
								v = max(abs(fitC(3))/h/tau(end),abs(vz));  vt = sqrt(v^2-vz^2);
							else
								vt = abs(fitC(3))/hx/tau(end);  v = sqrt(vt^2+vz^2);
							end
						end
						
						Ms(iz,ix,iy,ja) = ms;  Mf(iz,ix,iy,ja) = mf;  D(iz,ix,iy,ja) = d;  V(iz,ix,iy,ja) = v;
						if (v > 0)  A(iz,ix,iy,ja) = vz/v;  end;
						GGf(iz,ix,iy,2:end,ja) = ms + mf*exp(-hx^2*vt^2*t.^2-h^2*vz^2*t.^2-q^2*d*t).*exp(1i*q*vz*t);  GGf(iz,ix,iy,1,ja) = 1;
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
