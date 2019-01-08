% find Y = A*x^3 + B*x^2 + C*x + D
function [A B C D R2] = FitPolyn3(Y,x,dim)

ConvertMMA2;

if nargin < 3
	dim = 1;
end

	n = length(x);
	sz = size(Y);
	if sz(dim) ~= n
		error('sz(dim) ~= length(x)');
	end

	% make Y to [n m] (m = #channel)
	if dim > 1
		Y = shiftdim(Y,dim-1);
	end
	sz2 = size(Y);
	Y = Y(:,:);		% [n m]
	m = size(Y,2);

	% Y = X * C0
	X = [x(:).^3 x(:).^2 x(:) ones(n,1)];	% [n 4]
	C0 = inv(X'*X)*(X'*Y);			% [4 m]
	R2 = 1 - mean((Y - X*C0).^2,1) ./ mean((Y - repmat(mean(Y,1),[n 1])).^2,1);		% [1 m]

	A = C0(1,:);  B = C0(2,:);  C = C0(3,:);  D = C0(4,:);  R2 = R2(:);		% [m 1]

	% make A B C to size(Y) where size(Y,dim) = 1
	sz3 = sz2;  sz3(1) = 1;
	A = reshape(A,sz3);  B = reshape(B,sz3);  C = reshape(C,sz3);  D = reshape(D,sz3);  R2 = reshape(R2,sz3);
	if dim > 1
		A = shiftdim(A,length(sz)-dim+1);
		B = shiftdim(B,length(sz)-dim+1);
		C = shiftdim(C,length(sz)-dim+1);
		D = shiftdim(D,length(sz)-dim+1);
		R2 = shiftdim(R2,length(sz)-dim+1);
	end
