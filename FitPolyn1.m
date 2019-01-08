%% find Y = A*x + B
% bBzero : force B=0?

function [A,B,R2] = FitPolyn1(Y,x,dim,bBzero)

ConvertMMA2;

if nargin < 4
    bBzero = false;
end
if nargin < 3
	dim = 1;
end

% tic;

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
    if bBzero
        X = x(:);	% [n 1]
        C0 = inv(X'*X)*(X'*Y);			% [1 m]
        R2 = 1 - mean((Y - X*C0).^2,1) ./ mean((Y - repmat(mean(Y,1),[n 1])).^2,1);		% [1 m]

        A = C0(1,:);  B = zeros(1,m);  R2 = R2(:);		% [m 1]
    else
        X = [x(:) ones(n,1)];	% [n 2]
        C0 = inv(X'*X)*(X'*Y);			% [2 m]
        R2 = 1 - mean((Y - X*C0).^2,1) ./ mean((Y - repmat(mean(Y,1),[n 1])).^2,1);		% [1 m]

        A = C0(1,:);  B = C0(2,:);  R2 = R2(:);		% [m 1]
    end
        
	% make A B C to size(Y) where size(Y,dim) = 1
	sz3 = sz2;  sz3(1) = 1;
	A = reshape(A,sz3);  B = reshape(B,sz3);  R2 = reshape(R2,sz3);
	if dim > 1
		A = shiftdim(A,length(sz)-dim+1);
		B = shiftdim(B,length(sz)-dim+1);
		R2 = shiftdim(R2,length(sz)-dim+1);
    end
    
% disp([datestr(now,'HH:MM') '  FitPolyn1 (' num2str(toc) ' s)']);

