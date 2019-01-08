function ret = Coregister2D_xcr(A, A0, r1, cor1, t1, id)

            [nx ny] = size(A);
            A1 = RotateImage(A,r1,cor1,t1);
                    %{
                    A1 = RotateImage(A(:,200:end),pi/20,[300 0],[100 0]);
                    figure(1);  clf;
                        subplot(121);  hold on;  imagesc(A(:,200:end)');  axis image;
                        subplot(122);  hold on;  imagesc(A1');  axis image;
                    %}
            sqA = A(:);  sqA0 = A0(:);  sqA1 = A1(:);
            idx = find(~isnan(sqA1));
            M = zeros(nx,ny);  M = M(:);  M(idx) = 1;  M = reshape(M,[nx ny]);
            xcr0 = GetXcr(sqA(idx),sqA0(idx));
            xcr = GetXcr(sqA1(idx),sqA0(idx));
                if id > 0
                    cm = zeros(128,3);  cm(1:64,:) = gray(64);  cm(65:128,1) = linspace(0,1,64);
%                     figure(figno);  clf;  colormap(cm);  hold on; 
                    subplot(4,4,[1 2 5 6]);  cla;  colormap(gca,cm);  hold on; 
                        PlotImage(Rescale(A0',[0 0.9],[MeanLow(A0,0.1) MeanHigh(A0,0.9)]),0,[0 2],0);
                        img = Rescale(A1'.*M',[0.1 0.7],[MeanLow(sqA1(idx),0.1) MeanHigh(sqA1(idx),0.9)]);
                        PlotImage(img+1,0,[0 2],0,[],img);
                        set(gca,'CLim',[0 2]);
						title(['xcr improving from ' num2str(xcr0,3) ' to ' num2str(xcr,3)]);
                    pause(.1);
                end
    
ret = xcr;
       