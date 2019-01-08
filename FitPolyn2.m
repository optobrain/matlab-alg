% find Y = A*x^2 + B*x + C
%{
function [A B C E] = FitPolyn2(Y,x,dim)

	sz = size(Y);
	szX0 = ones(1,length(sz));  szX0(dim) = length(x);
	szX = sz;  szX(dim) = 1;
	X = repmat(reshape(x,szX0),szX);
		
	a = sum(X.^4,dim);  b = sum(X.^3,dim);  c = sum(X.^2,dim);
	d = sum(X.^3,dim);  e = sum(X.^2,dim);  f = sum(X,dim);
	g = sum(X.^2,dim);  h = sum(X,dim);  k = ones(size(h));

	deter = a.*(e.*k-f.*h) + b.*(f.*g-k.*d) + c.*(d.*h-e.*g);
	A = e.*k - f.*h;  D = c.*h - b.*k;  G = b.*f - c.*e;
	B = f.*g - d.*k;  E = a.*k - c.*g;  H = c.*d - a.*f;
	C = d.*h - e.*g;  F = g.*b - a.*h;  K = a.*e - b.*d;

%	[a b c; d e f; g h k]
%	inv([a b c; d e f; g h k])
%	[A D G; B E H; C F K]/deter
%	[A D G; B E H; C F K]/deter * [a b c; d e f; g h k]
    
	AA = ( A.*sum(Y.*X.^2,dim) + D.*sum(Y.*X,dim) + G.*sum(Y,dim) ) ./ deter;
	BB = ( B.*sum(Y.*X.^2,dim) + E.*sum(Y.*X,dim) + H.*sum(Y,dim) ) ./ deter;
	CC = ( C.*sum(Y.*X.^2,dim) + F.*sum(Y.*X,dim) + K.*sum(Y,dim) ) ./ deter;

	A = AA;  B = BB;  C = CC;

	E = sum(abs( repmat(A,szX0).*X.^2 + repmat(B,szX0).*X + repmat(C,szX0) - Y ).^2,dim);

%}


% find Y = A*x^2 + B*x
%{
function [A B E] = FitPolyn2(Y,x,dim)

	sz = size(Y);
	szX0 = ones(1,length(sz));  szX0(dim) = length(x);
	szX = sz;  szX(dim) = 1;
	X = repmat(reshape(x,szX0),szX);
		
	den = sum(X.^4,dim).*sum(X.^2,dim) - sum(X.^3,dim).^2;
	A = ( sum(X.^2,dim).*sum(X.^2.*Y,dim) - sum(X.^3,dim).*sum(X.*Y,dim) ) ./ den;
	B = ( sum(X.^4,dim).*sum(X.*Y,dim) - sum(X.^3,dim).*sum(X.^2.*Y,dim) ) ./ den;

	E = sum(abs( repmat(A,szX0).*X.^2 + repmat(B,szX0).*X - Y ).^2,dim);
%}

% find Y = A*x^2 + B*x + C
function [A B C R2] = FitPolyn2(Y,x,dim)

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
	X = [x(:).^2 x(:) ones(n,1)];	% [n 3]
	C0 = inv(X'*X)*(X'*Y);			% [3 m]
	R2 = 1 - mean((Y - X*C0).^2,1) ./ mean((Y - repmat(mean(Y,1),[n 1])).^2,1);		% [1 m]

	A = C0(1,:);  B = C0(2,:);  C = C0(3,:);  R2 = R2(:);		% [m 1]

	% make A B C to size(Y) where size(Y,dim) = 1
	sz3 = sz2;  sz3(1) = 1;
	A = reshape(A,sz3);  B = reshape(B,sz3);  C = reshape(C,sz3);  R2 = reshape(R2,sz3);
	if dim > 1
		A = shiftdim(A,length(sz)-dim+1);
		B = shiftdim(B,length(sz)-dim+1);
		C = shiftdim(C,length(sz)-dim+1);
		R2 = shiftdim(R2,length(sz)-dim+1);
	end
