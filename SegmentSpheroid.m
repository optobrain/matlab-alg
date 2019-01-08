
% II can be Mf, so I didn't use log10 here.  You have to use log10 in the main script

function [iz,ix,iy,rx,rz] = SegmentSpheroid(II,ix,iy,rx,rz)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jonghwan_lee@brown.edu');
end

    [nz,nx,ny] = size(II);
    
                    xx = max(ix-rx*2,1):min(ix+rx*2,nx);  yy = max(iy-rx*2,1):min(iy+rx*2,ny);
                    Iroi = II(:,xx,yy);
                    
                    % from MIP, quantify (error pixel number) = (I<Ithr inside the circle) + (I>Ithr outside the circle)
                    Isph = squeeze(max(Iroi,[],1));  Isq = Isph(:);  idx = kmeans(Isq,2);  thr = min(Isq(idx==2));  Ib = (Isph >= thr);  
                    if sum(Ib(:)) < numel(Ib)*0.1 || sum(Ib(:)) > numel(Ib)*0.9
%                         sum(Ib(:))/numel(Ib)
                        thr = GetSorted(Isq,0.75);  Ib = (Isph >= thr);
                    end
%                         subplot(4,5,4);  hold off;  cla;  hold on;  
%                             PlotImage(squeeze(max(Iroi,[],1))',false,limI);  xlabel('X');  ylabel('Y');  
%                             line(ix-xx(1)+1,iy-yy(1)+1,'color','r','marker','+');
%                             t = (-pi:0.1:pi);  x = ix-xx(1)+1+rx*cos(t);  y = iy-yy(1)+1+rx*sin(t);  line(x,y,'color','r');
%                         subplot(4,5,4);  hold off;  cla;  axis auto;  hist(Isq);  line(thr*[1 1],get(gca,'YLim'));  xlabel('I');
                        subplot(4,5,4);  hold off;  cla;  hold on;  
                            PlotImage(Ib',false,[0 1]);  xlabel('X');  ylabel('Y');  
                            line(ix-xx(1)+1,iy-yy(1)+1,'color','r','marker','+');
                            t = (-pi:0.1:pi);  x = ix-xx(1)+1+rx*cos(t);  y = iy-yy(1)+1+rx*sin(t);  line(x,y,'color','r');
                    [gy,gx] = meshgrid(1:length(yy),1:length(xx));
                    r = (gx-(ix-xx(1)+1)).^2/rx^2 + (gy-(iy-yy(1)+1)).^2/rx^2;  
                    err1 = sum(sum(r<=1))-sum(sum(Ib(r<=1))) + sum(sum(Ib(r>1)));  
                    
                    % from MIP, find the center (x,y) and rx by minimizing (error pixel number) 
                    iix = ix+(-4:4);  iiy = iy+(-4:4);  irx = rx+(-4:4);
                    err = zeros(length(iix),length(iiy),length(irx));
                    for i1=1:length(iix)
                    for i2=1:length(iiy)
                    for i3=1:length(irx)
                        r = (gx-(iix(i1)-xx(1)+1)).^2/irx(i3)^2 + (gy-(iiy(i2)-yy(1)+1)).^2/irx(i3)^2;  
                        err(i1,i2,i3) = sum(sum(r<=1))-sum(sum(Ib(r<=1))) + sum(sum(Ib(r>1)));
                    end
                    end
                    end
                    [m,im] = FindMax(1./err);
                    ix = iix(im(1));  iy = iiy(im(2));  rx = irx(im(3));
                            line(ix-xx(1)+1,iy-yy(1)+1,'color','b','marker','+');
                            t = (-pi:0.1:pi);  x = ix-xx(1)+1+rx*cos(t);  y = iy-yy(1)+1+rx*sin(t);  line(x,y,'color','b');
                            title(['error: ' num2str(err1) ' => ' num2str(1/m)]);
                    
                    % find the center z where I>Ithr is maximum
                    Isq = Iroi(:);  idx = kmeans(Isq,2);  thr = min(Isq(idx==2));  Ib = (Iroi >= thr);  
                    if sum(Ib(:)) < numel(Ib)*0.01 || sum(Ib(:)) > numel(Ib)*0.2
%                         sum(Ib(:))/numel(Ib)
                        thr = GetSorted(Isq,0.9);  Ib = (Iroi >= thr);
                    end
                    Iz = convn(mean(Iroi(:,:),2),ones(3,1)/3,'same');  Ibz = convn(mean(Ib(:,:),2),ones(3,1)/3,'same');
%                         subplot(4,5,13);  hold off;  cla;  axis auto;  hist(Isq);  line(thr*[1 1],get(gca,'YLim'));  xlabel('I');
                        subplot(4,5,9);  hold off;  cla;  axis auto;  plotyy(1:nz,Iz,1:nz,Ibz);  
                    [m,iz] = max(Ibz);                     
                    
                    % from the top half volume, quantify (error voxel number) = (I<Ithr inside the ellipsoid) + (I>Ithr outside the ellipsoid)
                    zz = 1:iz;  Ib = Ib(zz,:,:);  % top half
                    [gx,gz,gy] = meshgrid(1:length(xx),1:iz,1:length(yy));
                    r = (gz-iz).^2/rz^2 + (gx-(ix-xx(1)+1)).^2/rx^2 + (gy-(iy-yy(1)+1)).^2/rx^2;  
                    err1 = sum(sum(sum(r<=1)))-sum(sum(sum(Ib(r<=1)))) + sum(sum(sum(Ib(r>1))));  
                        subplot(4,5,14);  hold off;  cla;  axis auto;  
                            PlotImage(max(Ib,[],3),false,[0 1]);  xlabel('X');  ylabel('Z');
                            t = (-pi:0.1:0);  x = ix-xx(1)+1+rx*cos(t);  z = iz+rz*sin(t);  line(x,z,'color','r');
                            
                    % from the top half, find rz by minimizing (error voxel number) 
                    irz = rz+(-floor(rz/2):(iz-rz));
%                     irz = rz+(-10:10);
                    if length(irz) > 1
                        err = zeros(length(irz),1);
                        for i1=1:length(irz)
                            r = (gz-iz).^2/irz(i1)^2 + (gx-(ix-xx(1)+1)).^2/rx^2 + (gy-(iy-yy(1)+1)).^2/rx^2;  
                            err(i1) = sum(sum(sum(r<=1)))-sum(sum(sum(Ib(r<=1)))) + sum(sum(sum(Ib(r>1))));  
                        end
                        [m,im] = max(1./err);
                        rz = irz(im);
                    end
                            t = (-pi:0.1:0);  x = ix-xx(1)+1+rx*cos(t);  z = iz+rz*sin(t);  line(x,z,'color','b');
                            title(['error: ' num2str(err1) ' => ' num2str(1/m)]);
                            
                            
                            