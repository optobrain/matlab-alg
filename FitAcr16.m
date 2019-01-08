%% Fit DLS after -COR

function [Mf,D,V,A,R,GGf] = FitAcr16(tau, GG, xzr, Tol)

ConvertMMA2;

if (nargin < 4) Tol = 1e-3;  end;
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
        na = 2;
		Mf = zeros(nz,nx,ny,na);  D = Mf;  V = Mf;  A = Mf;  R = Mf;  GGf = complex(ones(nz,nx,ny,ntau,na));
		t = reshape(tau(2:end),[ntau-1 1]);		
		tn = t / tau(end);						% tn(end) = 1 for TolX
        
		% guess the initial Mf0 with three points
		X = [tau(2)^2 tau(2) 1; tau(3)^2 tau(3) 1; tau(4)^2 tau(4) 1];
		C = Multiply(inv(X),shiftdim(abs(GG(:,:,:,2:4)),3),1);
		Mf0 = min(max( shiftdim(C(3,:,:,:),1) ,0),1);
		clear X C;

	
	%% Fit 1: Vz -> Mf & D & Vt
        ja = 1;
        Vz = V;

		% Vz
        bw = true;
        GGa = unwrap(angle(GG),[],4);
		T = repmat(reshape(tau,[1 1 1 ntau]),[nz nx ny 1]);
		if bw
			Vz = sum(GGa.*T.*abs(GG),4) ./ sum(T.^2.*abs(GG),4) /q;
		else
			Vz = sum(GGa.*T,4) ./ sum(T.^2,4) /q;		% = Multiply(GGa,tau') / (tau*tau') /q
		end
		clear GGa T;
        
        % Mf D Vt
		for iz=1:nz
			for ix=1:nx
				for iy=1:ny

					% prepare
					g = squeeze(GG(iz,ix,iy,2:end));
					di0 = .5/q^2/tau(end);
					vi0 = .5/h/tau(end);
                    vz = Vz(iz,ix,iy,ja);
                    vt0 = sqrt(max(vi0^2-vz^2,0));
					mf0 = Mf0(iz,ix,iy);
                    
                    % fit
                    fitE2 = @(c) sum( abs( min(abs(c(1)),1).*exp( -h^2*vz^2*t.^2 -c(2).^2.*tn.^2 -abs(c(3)).*tn ) - abs(g) ).^2 );
                    fitC0 = [mf0 hx*vt0*tau(end) q^2*di0*tau(end)];
                    fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolFun',Tol,'TolX',Tol));
                    mf = min(abs(fitC(1)),1);  vt = abs(fitC(2))/hx/tau(end);  d = abs(fitC(3))/q^2/tau(end);
                    vt = abs(fitC(3))/hx/tau(end);  v = sqrt(vt^2+vz^2);

                    % 
                    Mf(iz,ix,iy,ja) = mf;  D(iz,ix,iy,ja) = d;  V(iz,ix,iy,ja) = v;
                    if (v > 0)  A(iz,ix,iy,ja) = vz/v;  end;
                    GGf(iz,ix,iy,2:end,ja) = mf*exp(-hx^2*vt^2*t.^2-h^2*vz^2*t.^2-q^2*d*t).*exp(1i*q*vz*t);  GGf(iz,ix,iy,1,ja) = 1;
                    R(iz,ix,iy,ja) = 1 - sum( abs(GG(iz,ix,iy,:)-GGf(iz,ix,iy,:,ja)).^2 ,4) / sum( abs(GG(iz,ix,iy,:)-mean(GG(iz,ix,iy,:),4)).^2 ,4);									
                    
                end
            end
        end
        
                    
	%% Fit 2: Mf => D & Vt & Vz
        ja = 2;
        bw = true;
        
        % Mf => D Vt Vz
		for iz=1:nz
			for ix=1:nx
				for iy=1:ny

					% prepare
					g = squeeze(GG(iz,ix,iy,2:end));
					mf0 = Mf0(iz,ix,iy);

                    % fit: Mf first and then the others
                    fitE2 = @(c) sum( abs( min(abs(c(1)),1) .* FitAcr14exp(tn, g, 0, min(abs(c(1)),1), q,h,bw) - g ).^2 );
                    fitC0 = [mf0];
                    fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolFun',Tol,'TolX',Tol));
                    mf = min(abs(fitC(1)),1);  
                    [~,d,vt,vz] = FitAcr14exp(tn, g, 0, mf, q,h,bw);
                    d = d/q^2/tau(end);  vt = vt/hx/tau(end);  vz = vz/q/tau(end);  v = sqrt(vt^2+vz^2);

                    Mf(iz,ix,iy,ja) = mf;  D(iz,ix,iy,ja) = d;  V(iz,ix,iy,ja) = v;
                    if (v > 0)  A(iz,ix,iy,ja) = vz/v;  end;
                    GGf(iz,ix,iy,2:end,ja) = mf*exp(-hx^2*vt^2*t.^2-h^2*vz^2*t.^2-q^2*d*t).*exp(1i*q*vz*t);  GGf(iz,ix,iy,1,ja) = 1;
                    R(iz,ix,iy,ja) = 1 - sum( abs(GG(iz,ix,iy,:)-GGf(iz,ix,iy,:,ja)).^2 ,4) / sum( abs(GG(iz,ix,iy,:)-mean(GG(iz,ix,iy,:),4)).^2 ,4);									
                    
                end
            end
        end
    
        
                          