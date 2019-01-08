%% Get PCA with pca()

% Y [nch nt]
% Wa : for inverse

function [Yc, W] = GetPCA(Y, Wa)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jlee@optobrain.com');
end

if nargin < 2
    W = pca(Y');
    Yc = (Y'*W)';
    
else  % inverse
    W = Wa;
    Yc = (Y'*inv(W))';
end
