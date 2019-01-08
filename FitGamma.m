%% Fit multichannel data [nt nch] with my custom Gamma function
% Y = A (t-T)^m exp(-B(t-T)) * step(t-T),  B > 0
% Peak & peak time
%		m = 1: tp = T + 1/A/B,  Yp = 1/B*exp(-1/A)
%		m > 1: tp = T + m/B,  Yp = A(m/B)^m exp(-m)

% tt : time array
% Y [nt nch]
% m : order of Gamma function
% cA : array of initial guess of A
% cB : array of initial guess of B
% cT : array of initial guess of T
% bTpos : T is always > 0 ?
% Bmax : set max of B to prevent very sharp fitting

% A, B, T, R, P : fitted coefficients and p-value with max R among candidate initial guesses
% tp, Yp : peak time and amplitudes

function [A,B,T,R,P,tp,Yp,Yfit] = FitGamma(tt, Y, m, cA, cB, cT, bTpos, Bmax)

ConvertMMA2;

if nargin < 8
	Bmax = inf;
end
if nargin < 7
	bTpos = 0;
end

	[nt,nch] = size(Y);
	nA = length(cA);  nB = length(cB);  nT = length(cT);

	A = zeros(1,nch);  B = A;  T = A;  R = A;  P = A;
	for ich=1:nch
		y = Y(:,ich);
		a = zeros(nA,nB,nT);  b = a;  t = a;  r = a;
		for ia=1:nA
			for ib=1:nB
				for it=1:nT
					if bTpos == 1		% T > 0
						fitE2 = @(c) sum( ( c(1) .* (tt-abs(c(3))).^m .* exp(-min(abs(c(2)),Bmax).*(tt-abs(c(3)))) .* (1+sign(tt-abs(c(3))))/2  - y ).^2 );
						fitC0 = [cA(ia) cB(ib) cT(it)];
						fitC = fminsearch(fitE2, fitC0, optimset('Display','off'));
						a(ia,ib,it) = fitC(1);  b(ia,ib,it) = min(abs(fitC(2)),Bmax);  t(ia,ib,it) = abs(fitC(3));
					else
						fitE2 = @(c) sum( ( c(1) .* (tt-c(3)).^m .* exp(-min(abs(c(2)),Bmax).*(tt-c(3))) .* (1+sign(tt-c(3)))/2  - y ).^2 );
						fitC0 = [cA(ia) cB(ib) cT(it)];
						fitC = fminsearch(fitE2, fitC0, optimset('Display','off'));
						a(ia,ib,it) = fitC(1);  b(ia,ib,it) = min(abs(fitC(2)),Bmax);  t(ia,ib,it) = fitC(3);
					end
					yf = a(ia,ib,it)*(tt-t(ia,ib,it)).^m.*exp(-b(ia,ib,it)*(tt-t(ia,ib,it))).*(1+sign(tt-t(ia,ib,it)))/2;
					r(ia,ib,it) = 1 - mean(( yf -y).^2) / mean(( mean(y) -y).^2);
				end
			end
		end
		[r,im] = FindMax(r);
        if length(im) == 3
    		ia = im(1);  ib = im(2);  it = im(3);
        elseif length(im) == 2
    		ia = im(1);  ib = im(2);  it = 1;
        else
    		ia = im(1);  ib = 1;  it = 1;
        end            
		a = a(ia,ib,it);  b = b(ia,ib,it);  t = t(ia,ib,it);
		yf = a*(tt-t).^m.*exp(-b*(tt-t)).*(1+sign(tt-t))/2;
		[be,bint,r2,rint,stat] = regress(y(:),[ones(size(yf(:))) yf(:)]);  
		p = stat(3);

		A(ich) = a;  B(ich) = b;  T(ich) = t;  R(ich) = r;  P(ich) = p;  Yfit(:,ich) = yf;

		if mod(ich,100) == 0
			disp(['... FitGamma ' num2str(ich) '/' num2str(nch) ' ' datestr(now,'HH:MM')]);
		end
	end

%%		m = 1: tp = T + 1/A/B,  Yp = 1/B*exp(-1/A)
%%		m > 1: tp = T + m/B,  Yp = A(m/B)^m exp(-m)

	if m == 1
		tp = T + 1./A./B; 
		Yp = 1./B.*exp(-1./A);
	else
		tp = T + m./B;
		Yp = A.*(m./B).^m.*exp(-m);
	end






