%% Fit DLS only for D (V=0)
% xzr = X_resol / Z_resol where Z_resol = 6.5 um AIR (FWHM)
% bms0 = force Ms=0?  This was more accurate.
% balp = force alpha = 1 ?  
% blsq = use lsqcurvefit instead of fminsearch ? 

function [Ms,Mf,D,ALP,R,GGf] = FitAcr17(tau, GG, xzr, bms0, balp, blsq, Tol)

ConvertMMA2;

if (nargin < 7) Tol = 1e-3;  end;
if (nargin < 6) blsq = true;  end;  % lsqcurvefit resulted in slightly higher R2
if (nargin < 5) balp = true;  end;  % balp=false is very unstable (R2 is very low)
if (nargin < 4) bms0 = true;  end;
if (nargin < 3) xzr = 1;  end;

	%% Prepare

		% constant
        % for simulation
		k0 = 2*pi/1.31e-6;
%         k_range = 0.17e-6*2.2;  % ideal: 2.2 * FWHM
% 		dk = ( 2*pi/(1.31e-6-0.17e-6/2)-2*pi/(1.31e-6+0.17e-6/2) )/2*sqrt(2*log(2));
		n = 1.35;  % n=1.35 resulted in 2x higher D than true values
		q = 2*n*k0;
% 		dq = 2*n*dk;
% 		h = dq/2;
% 		hx = h/xzr;
        
		% prepare
		[nz,nx,ny,ntau] = size(GG);
% 		dt = mean(diff(tau));
        na = 1;
		Ms = zeros(nz,nx,ny,na);  Mf = Ms;  D = Mf;  ALP = Mf;  R = Mf;  GGf = ones(nz,nx,ny,ntau,na);
		t = reshape(tau(2:end),[ntau-1 1]);		
		tn = t / tau(end);						% tn(end) = 1 for TolX
        
		% guess the initial Mf0 with three points
		X = [tau(2)^2 tau(2) 1; tau(3)^2 tau(3) 1; tau(4)^2 tau(4) 1];
		C = Multiply(inv(X),shiftdim(abs(GG(:,:,:,2:4)),3),1);
		Me0 = 1 - min(max( shiftdim(C(3,:,:,:),1) ,0),1);
		clear X C;
        Mf0 = (1-Me0)/2;  Ms0 = Mf0;

	
	%% Fit 1: Ms & Mf & D
        ja = 1;
        
        % Mf D Vt
		for iz=1:nz
			for ix=1:nx
				for iy=1:ny

					% prepare
					g = squeeze(abs(GG(iz,ix,iy,2:end)));
                    ms0 = Ms0(iz,ix,iy);
					mf0 = Mf0(iz,ix,iy);
                    
                    % fit 
                    if balp  % force alpha=1
                        if bms0  % force Ms=0, which was more accurate
                            if blsq
                                fitF = @(c,tn) c(1).*exp( -c(2).*tn );
                                fitC0 = [mf0*2 0.5];  % di0 = 0.5/q^2/tau(end) corresponds to 10 um2/s
                                fitC = lsqcurvefit(fitF,fitC0,tn,g,[0 0],[1 5], optimset('Display','off','TolFun',Tol,'TolX',Tol));
                                ms = 0;  mf = fitC(1);  d = fitC(2)/2/q^2/tau(end);  
                            else
                                fitE2 = @(c) sum( ( min(abs(c(1)),1).*exp( -abs(c(2)).*tn ) - g ).^2 );
                                fitC0 = [mf0*2 0.5];
                                fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolFun',Tol,'TolX',Tol));
                                ms = 0;  mf = min(abs(fitC(1)),1);  d = abs(fitC(2))/2/q^2/tau(end);
                            end
                        else
                            fitE2 = @(c) sum( ( min(abs(c(1)),1) + min(abs(c(2)),1-min(abs(c(1)),1)).*exp( -abs(c(3)).*tn ) - g ).^2 );
                            fitC0 = [ms0 mf0 0.5];
                            fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolFun',Tol,'TolX',Tol));
                            ms = min(abs(fitC(1)),1);  mf = min(abs(fitC(2)),1-ms);  d = abs(fitC(3))/2/q^2/tau(end);
                        end
                        alp = 1;
                    else
                        if bms0  % force Ms=0, which was more accurate
                            if blsq
                                fitF = @(c,tn) c(1).*exp( -c(2).*tn.^c(3) );
                                fitC0 = [mf0*2 0.5 1];
                                fitC = lsqcurvefit(fitF,fitC0,tn,g,[0 0 0],[1 5 2], optimset('Display','off','TolFun',Tol,'TolX',Tol));
                                ms = 0;  mf = fitC(1);  d = fitC(2)/2/q^2/tau(end);  alp = fitC(3);

%                                 fitF = @(c,tn) c(1) -c(2).*tn.^c(3);
%                                 fitC0 = [log(mf0*2) 0.5 1.5]
%                                 fitC = lsqcurvefit(fitF,fitC0,tn,log(g),[log(0.1) 0 0],[0 5 2])
%                                 ms = 0;  mf = exp(fitC(1));  d = fitC(2)/2/q^2/tau(end);  alp = fitC(3);
% 
%                                 fitF = @(c,tn) -c(1).*tn.^c(2);
%                                 fitC0 = [0.5 1]
%                                 fitC = lsqcurvefit(fitF,fitC0,tn,log(g),[0 0],[5 2])
%                                 ms = 0;  mf = 1;  d = fitC(1)/2/q^2/tau(end);  alp = fitC(2);
% 
%                                 fitF = @(c,tn) c(1) + c(2).*log(tn);
%                                 fitC0 = [log(0.5) 1]
%                                 fitC = lsqcurvefit(fitF,fitC0,tn,log(-log(g)),[log(0.5*0.1) 0],[log(5) 2])
%                                 ms = 0;  mf = 1;  d = exp(fitC(1)/2/q^2/tau(end));  alp = fitC(2);
                            else
                                fitE2 = @(c) sum( ( min(abs(c(1)),1).*exp( -abs(c(2)).*tn.^min(abs(c(3)),2) ) - g ).^2 );
                                fitC0 = [mf0*2 0.5 1];
                                fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolFun',Tol,'TolX',Tol));
                                ms = 0;  mf = min(abs(fitC(1)),1);  d = abs(fitC(2))/2/q^2/tau(end);  alp = min(abs(fitC(3)),2);
                            end

                        else
                            fitE2 = @(c) sum( ( min(abs(c(1)),1) + min(abs(c(2)),1-min(abs(c(1)),1)).*exp( -abs(c(3)).*tn.^min(abs(c(4)),2) ) - g ).^2 );
                            fitC0 = [ms0 mf0 0.5 1];
                            fitC = fminsearch(fitE2, fitC0, optimset('Display','off','TolFun',Tol,'TolX',Tol));
                            ms = min(abs(fitC(1)),1);  mf = min(abs(fitC(2)),1-ms);  d = abs(fitC(3))/2/q^2/tau(end);  alp = min(abs(fitC(4)),2);
                        end
                    end                      

                    Ms(iz,ix,iy,ja) = ms;  Mf(iz,ix,iy,ja) = mf;  D(iz,ix,iy,ja) = d;  ALP(iz,ix,iy,ja) = alp;
                    GGf(iz,ix,iy,2:end,ja) = ms + mf*exp(-q^2*d*t.^alp);  GGf(iz,ix,iy,1,ja) = 1;
                    R(iz,ix,iy,ja) = 1 - sum( abs(abs(GG(iz,ix,iy,:))-GGf(iz,ix,iy,:,ja)).^2 ,4) / sum( abs(abs(GG(iz,ix,iy,:))-mean(abs(GG(iz,ix,iy,:)),4)).^2 ,4);									
                    
                end
            end
        end
        
                    
	%% Fit 2: Mf => D & Vt & Vz
    %{
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
    %}
    
        
                          