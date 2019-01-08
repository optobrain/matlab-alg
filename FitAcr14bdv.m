%% for ver. 14b :: fixing Vt0
% bw : weight of exponential

function [D,R] = FitAcr14bdv(tau, GG, Ms, Mf, q, h, hx, vt0)
if (nargin < 8)  vt0 = 0;   end;

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

	% Vz, Vt
        Vz = zeros(nz,nx,ny,1,nm);
        Vt = ones(nz,nx,ny,1,nm)*vt0;

	% D
		Gf = abs(G-Ms) ./ Mf .* exp(hx^2*repmat(Vt.^2,[1 1 1 nt 1]).*T.^2);  % = exp(-q^2*D*T)
        D = reshape(FitPolyn1(log(Gf),q^2*t,4,true),[nz nx ny 1 nm]);

	% R
		GGf = ones(nz,nx,ny,nt+1,nm);
		GGf(:,:,:,2:end,:) = Ms + Mf .* exp( -h^2*repmat(Vz.^2,[1 1 1 nt 1]).*(T*tau(end)).^2 -hx^2*repmat(Vt.^2,[1 1 1 nt 1]).*(T*tau(end)).^2 -q^2*repmat(D,[1 1 1 nt 1]).*T*tau(end) ) .* exp( 1i*q*repmat(Vz,[1 1 1 nt]).*T*tau(end));
		R = 1 - sum(abs(repmat(GG,[1 1 1 1 nm])-GGf).^2,4) ./ repmat(sum(abs(GG-repmat(mean(GG,4),[1 1 1 nt+1])).^2,4),[1 1 1 1 nm]);



