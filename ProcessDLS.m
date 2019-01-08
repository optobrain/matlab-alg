% thrMf: Mf<thrMf is converted to zero (default 0)
% thrR: voxels with R2<thrR will be ignored in medfilt (default 0)
% medfiltsize: size of the edge of the median filter kernel (default 3) :: should be an odd number
% GG and GGf can be [] when want to minimize the memory usage, and output can be [Ms,Mf,D,V,A,R,~,~]

function [Ms,Mf,D,V,A,R,GG,GGf] = ProcessDLS(Ms0,Mf0,D0,V0,A0,R0,GG0,GGf0,thrMf,thrR,medfiltsize)

ConvertMMA2;

if (nargin < 11)  medfiltsize = 3;  end;
if (nargin < 10)  thrR = 0;  end;
if (nargin < 9)  thrMf = 0;  end;

    %% Mf
    if thrMf > 0
        Mf0(Mf0<thrMf) = 0;
    end
    MfD0 = Mf0.*D0;
    
    [nz,nx,ny] = size(Ms0);
    Ms = Ms0;  Mf = Mf0;
    D = D0;  V = V0;  A = A0;  R = R0;
    GG = GG0;  GGf = GGf0;
    
    %% medfilt & R2
    nm = medfiltsize;  ne = (nm-1)/2;  im = (nm^3+1)/2;
    if ne > 0
        [x,z,y] = meshgrid(1:nm,1:nm,1:nm);
        x = x(:);  z = z(:);  y = y(:);
        parfor iz=(1+ne):nz-ne
            for iy=(1+ne):ny-ne
                for ix=(1+ne):nx-ne
                    zz = iz-ne:iz+ne;  xx = ix-ne:ix+ne;  yy = iy-ne:iy+ne;
                    mfd = MfD0(zz,xx,yy);
                    [~,is] = sort(mfd(:));  jz = zz(z(is(im)));  jx = xx(x(is(im)));  jy = yy(y(is(im)));
                    if thrR > 0
                        r = R0(zz,xx,yy);  mfd = mfd(:);  r = r(:);  ir = find(r>thrR);  
                        if isempty(ir)  % when no voxel has r>0.95, then only median without considering R2
                            ir = [1:nm^3];  
                        end
                        nr = length(ir);  jr = round((nr+1)/2);
                        [~,is] = sort(mfd(ir));  jz = zz(z(is(jr)));  jx = xx(x(is(jr)));  jy = yy(y(is(jr)));
                    end
%                     MfD(iz,ix,iy) = MfD0(jz,jx,jy);
                    Ms(iz,ix,iy) = Ms0(jz,jx,jy);  Mf(iz,ix,iy) = Mf0(jz,jx,jy);
                    D(iz,ix,iy) = D0(jz,jx,jy);  V(iz,ix,iy) = V0(jz,jx,jy);  A(iz,ix,iy) = A0(jz,jx,jy);  R(iz,ix,iy) = R0(jz,jx,jy);  
                    if numel(GG0) > 0
                        GG(iz,ix,iy,:) = GG0(jz,jx,jy,:);  
                    end
                    if numel(GGf0) > 0
                        GGf(iz,ix,iy,:) = GGf0(jz,jx,jy,:);  
                    end
                end
            end
        end
    end
    
        
    

