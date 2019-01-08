%% get cross section of DD at r=[iz ix iy] over the plane normal to the vector ev with the size of (2*ned+1) x (2*ned+1)
% need to assign figure window or subplot before calling this function

function ret = GetCrossSection(DD, r, ev, ned, bView)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jonghwan_lee@brown.edu');
end

if nargin < 5
    bView = false;
end

    iz = r(1);  ix = r(2);  iy = r(3);
    
    nb = 2*ned+1;
    [xb,yb,zb] = meshgrid(-ned:ned,-ned:ned,-ned:ned);	

    Db = shiftdim(DD(iz+(-ned:ned),ix+(-ned:ned),iy+(-ned:ned)),1);	 Db = double(Db);% [x y z]
    if ev(2) ==0 && ev(3) == 0  % the below results in NaN when ev is only nonzero in Z
        imagesc(Db(:,:,ned+1));  axis image;
        ret = Db(:,:,ned+1);
    else
        hsp = surf(-ned:ned,-ned:ned,zeros(nb));
        theta = acos(ev(1))*180/pi;  phi = atan(ev(3)/(ev(2)+eps))*180/pi;
        rotate(hsp,[0 0 1],-phi);  rotate(hsp,[-ev(2) ev(3) 0],theta);
        xd = get(hsp,'XData');  yd = get(hsp,'YData');  zd = get(hsp,'ZData');  
        hs = slice(xb,yb,zb,Db,xd,yd,zd);  set(gca,'ZDir','reverse');  axis image;  
        if bView
            view(-phi,90-theta);	
        end
        ret = flipud(get(hs,'CData'));  
    end
