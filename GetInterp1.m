%% compare different methods of interp1 and choose the best (smallest error)

% Y = [nt nch]

function [Yi,bestm] = GetInterp1(t,Y,ti,bFig)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jlee@optobrain.com');
end

if nargin < 4
    bFig = 0;
end

    [nt nch] = size(Y);  nti = length(ti);
    
    sm = {'spline' 'pchip' 'linear'};
    
    YY = zeros(nti,nch,3);
    for is=1:3
        if ti(1) < t(1) || ti(end) > t(end)
            YY(:,:,is) = interp1(t,Y,ti,sm{is},'extrap');
        else
            YY(:,:,is) = interp1(t,Y,ti,sm{is});
        end
    end
    
    ES = zeros(nch,3);
    if nti > nt
        iti = zeros(1,nt);
        for ii=1:nt
            [m iti(ii)] = min(abs(ti-t(ii)));
        end
        ES = squeeze(mean((YY(iti,:,:) - repmat(Y, [1 1 3])).^2,1));
    else
        iti = zeros(1,nti);
        for ii=1:nti
            [m iti(ii)] = min(abs(t-ti(ii)));
        end
        ES = squeeze(mean((YY - repmat(Y(iti,:),[1 1 3])).^2,1));
    end

    Yi = zeros(nti,nch);
    for ich=1:nch
        [m im] = min(ES(ich,:));
        Yi(:,ich) = YY(:,ich,im);
        bestm{ich} = sm{im};
                if bFig
                    cl = lines(4);
                    figure(11);  clf;
                        line(t,Y(:,ich),'color',cl(1,:),'marker','o');
                        for is=1:3
                            line(ti,YY(:,ich,is),'color',cl(is+1,:));
                        end
                        legend({'raw' 'spline' 'pchip' 'linear'});
                        title(mat2str(ES(ich,:),2));
                    pause;
                end
    end
    
       
    
    