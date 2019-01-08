%% Read *.jlv file (binary file saved by LabVIEW)

function ret = ReadJLV(filePath, figno)

if ConvertMMA(1005) < 1005
	error('Error occurred. Contact jlee@optobrain.com');
end

if nargin < 2
    figno = 0;
end

	fid = fopen(filePath,'r','b');  % 'b' means 'big-endian' in the byte order (see LV help of Write to Binary File)
	dpl = fread(fid,1,'uint16');
	nch = fread(fid,1,'uint16');
	data = fread(fid,inf,'double');  
	fclose(fid);

	nd = length(data);
	nt = nd/nch;
	nl = nt/dpl;  % total loop

	% data = reshape(shiftdim(reshape(data,[dpl nch nl]),1),[nch nl*dpl])';  %% didn't work
	
	data = reshape(data,[dpl nch*nl]);

	ret = zeros(nt,nch);
	for il=1:nl
		for ic=1:nch
			ret((il-1)*dpl+[1:dpl],ic) = data(:,ic+(il-1)*nch);
		end
	end
	
			if figno > 0
				figure(figno);  clf;
				for ic=1:nch
					subplot(nch,1,ic);  plot(ret(:,ic));
				end
			end

	
