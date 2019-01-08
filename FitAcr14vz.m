function [Exp Vz] = FitAcr14vz(tau, GG, ms, mf, q)

	% repmat
		G = GG(:,:,:,2:end);
		[nz nx ny nt] = size(G);
		t = tau(2:end)/tau(end);
		T = repmat(reshape(t,[1 1 1 nt]),[nz nx ny 1]);

	% vz
		Ga = unwrap(angle(G-ms),[],4);
		if bw == 1
			Vz = 1/q * sum(Ga.*T.*abs(G-ms).^2,4) ./ sum(T.^2.*abs(G-ms).^2);
		else
			Vz = 1/q * sum(Ga.*T,4) ./ sum(t.^2);
		end

	% Exp
		Exp = exp(1i*q*repmat(Vz,[1 1 1 nt]).*T);
		Vz = Vz /tau(end);

