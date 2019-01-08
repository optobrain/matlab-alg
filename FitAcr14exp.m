
% optimized for the use in fminsearch
% bw : weight of exponential

function [Exp D Vt Vz] = FitAcr14exp(tn, g, ms, mf, q, h, bw)
if (nargin < 7)  bw = true;   end;

    g = g(:);  tn = tn(:);
	% vz
		ga = unwrap(angle(g-ms));
		if bw
			Vz = sum(ga.*tn.*abs(g-ms).^2) ./ sum(tn.^2.*abs(g-ms).^2);
		else
			Vz = sum(ga.*tn) ./ sum(tn.^2);
		end

	% d, v
		gf = abs(g-ms) / max(mf,eps) .* exp((h/q)^2*Vz^2*tn.^2);
		if bw
			a = sum(tn.^4.*gf.^2);  b = sum(tn.^3.*gf.^2);  d = sum(tn.^2.*gf.^2);
			e = sum(tn.^2.*gf.^2.*log(gf));  f = sum(tn.*gf.^2.*log(gf));
		else
			a = sum(tn.^4);  b = sum(tn.^3);  d = sum(tn.^2);
			e = sum(tn.^2.*log(gf));  f = sum(tn.*log(gf));
		end
		Det = a.*d-b.^2;  Det = max(abs(Det),eps).*sign(Det);
		Vt = 1 ./ Det .* ( +d.*e -b.*f );
		D = 1 ./ Det .* ( -b.*e +a.*f );

	% Exp
		Vt = sqrt(max(-Vt,0));
		D = max(-D,0);
		Exp = exp(-Vt^2*tn.^2 -D*tn) .* exp(-(h/q)^2*Vz^2*tn.^2) .* exp(1i*Vz*tn);


