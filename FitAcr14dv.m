
% bw : weight of exponential

function [D Vt Vz R] = FitAcr14dv(tau, GG, Ms, Mf, q, h, hx, bw)
if (nargin < 8)  bw = true;   end;

if length(Ms) ~= length(Mf)
	error('The lengths of Ms and Mf do not agree.');
end

	% repmat
		nm = length(Ms);
		G = GG(:,:,:,2:end);
		[nz nx ny nt] = size(G);
		G = repmat(G,[1 1 1 1 nm]);
		t = tau(2:end)/tau(end);
		T = repmat(reshape(t,[1 1 1 nt 1]),[nz nx ny 1 nm]);
		Ms = repmat(reshape(Ms,[1 1 1 1 nm]),[nz nx ny nt 1]);
		Mf = repmat(reshape(Mf,[1 1 1 1 nm]),[nz nx ny nt 1]);

	% vz
		Ga = unwrap(angle(G-Ms),[],4);
		if bw
			Vz = 1/q * sum(Ga.*T.*abs(G-Ms).^2,4) ./ sum(T.^2.*abs(G-Ms).^2,4);
		else
			Vz = 1/q * sum(Ga.*T,4) ./ sum(t.^2);
		end

	% D, V (normalized by tau(end))
		Gf = abs(G-Ms) ./ Mf .* exp(h^2*repmat(Vz.^2,[1 1 1 nt 1]).*T.^2);
		if bw
			a = sum(T.^4.*Gf.^2,4);  b = sum(T.^3.*Gf.^2,4);  d = sum(T.^2.*Gf.^2,4);
			e = sum(T.^2.*Gf.^2.*log(Gf),4);  f = sum(T.*Gf.^2.*log(Gf),4);
		else
			a = sum(T.^4,4);  b = sum(T.^3,4);  d = sum(T.^2,4);
			e = sum(T.^2.*log(Gf),4);  f = sum(T.*log(Gf),4);
		end
		Det = a.*d-b.^2;  Det = max(abs(Det),eps).*sign(Det);
		Vt = 1 ./ Det .* ( +d.*e -b.*f );
		D = 1 ./ Det .* ( -b.*e +a.*f );

		Vz = Vz /tau(end);
		Vt = sqrt(max(-Vt,0))/hx /tau(end);
		D = max(-D,0)/q^2 /tau(end);

	% R
		GGf = ones(nz,nx,ny,nt+1,nm);
		GGf(:,:,:,2:end,:) = Ms + Mf .* exp( -h^2*repmat(Vz.^2,[1 1 1 nt]).*(T*tau(end)).^2 -hx^2*repmat(Vt.^2,[1 1 1 nt]).*(T*tau(end)).^2 -q^2*repmat(D,[1 1 1 nt]).*T*tau(end) ) .* exp( 1i*q*repmat(Vz,[1 1 1 nt]).*T*tau(end));
		R = 1 - sum(abs(repmat(GG,[1 1 1 1 nm])-GGf).^2,4) ./ repmat(sum(abs(GG-repmat(mean(GG,4),[1 1 1 nt+1])).^2,4),[1 1 1 1 nm]);



