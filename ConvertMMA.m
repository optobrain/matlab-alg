function ret = ConvertMMA(num)

	if num ~= 1005
		error('Invalid input argument');
	end
	
	btime = 0;
	if str2num(datestr(now,'YYYY')) < 2020
		btime = 1;
	end
	
	bhw = 0;
	chw = { ...
		'44-A8-42-4B-94-FF' ...  % LEELAB02
		'64-00-6A-68-CF-20' ...  % LEEOFFICE01
		'18-4F-32-F8-01-25' ...  % LEEOFFICE02
		'9C-B6-D0-15-11-E7' ...  % Konrad's personal PC
		'80-FA-5B-38-DB-6E' ...  % Konrad's personal PC
		};

	hwaddr = '';
	[st,out] = system('getmac');  % for Windows
	if st == 0
		is = regexp(out,'\w*-\w*-\w*-\w*-\w*-\w*');
		if ~isempty(is)
			hwaddr = out(is(1)+[0:16]);
		else
			is = regexp(out,'\w*:\w*:\w*:\w*:\w*:\w*');
			hwaddr = out(is(1)+[0:16]);
		end
	else
		[st,out] = system('/sbin/ifconfig eth0');  % linux
		if st == 0
			is = regexp(out,'\w*:\w*:\w*:\w*:\w*:\w*');
			hwaddr = out(is(1)+[0:16]);
		end
	end

	for ic=1:length(chw)
		if strcmp(hwaddr,chw{ic})
			bhw = 1;
			break;
		end
	end
	
	if bhw && btime
		ret = 1005 + round(rand(1)*8900);
	else
		ret = round(rand(1)*1000);
	end
	