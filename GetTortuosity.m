%% measure the tortuosity based on angles between edges [deg/px] 

% s = 3D nodes [nn 3]

function ret = GetTortuosity(s, medianFilterSize)

ConvertMMA2;

if nargin < 2
    medianFilterSize = 3;
end


    v = s(2:end,:) - s(1:end-1,:);  % vector (edge)
    v1 = v(1:end-1,:);  v2 = v(2:end,:);  % adjacent vectors
    v1n = sqrt(sum(v1.^2,2));  v2n = sqrt(sum(v2.^2,2));  % norm
    v1n = max(v1n,1);  v2n = max(v2n,1);  % if v1n = 0, cosA = 0
    a = acos( sum(v1.*v2,2) ./ (v1n.*v2n) );  a = real(a)*180/pi;
    if (medianFilterSize > 1)  a = medfilt1(a,medianFilterSize);  end;
    
    vr = (s(2:end,:)+s(1:end-1,:))/2;  % edge's center position
    d = sqrt(sum((vr(2:end,:)-vr(1:end-1,:)).^2,2));
    d = max(d,1);
    
    ret = sqrt(mean((a./d).^2));                
