% step 2, read out data from eye tracker
clear all
format compact
fs = filesep;
dir_data = 'PATH';
addpath(genpath([dir_data 'Exp1']))
addpath(genpath([dir_data 'Exp2']))
addpath(genpath([dir_data fs 'matlab']));
addpath(genpath([dir_data fs 'matlab' fs 'analysis']));
cd(dir_data)

% Exp1
subj = { 'p1' 'p2'}

% Exp2

subj = {'ea01b' 'ea02b' 'ea03b' 'ea04b' 'ea05b' 'ea06b' 'ea07b' 'ea08b' 'ea09b' 'ea10b' ...
        'ea11b' 'ea12b' 'ea13b' 'ea14b' 'ea15b' 'ea16b' 'ea17b' 'ea18b' 'ea19b' 'ea20b' ...
        'ea21b' 'ea22b' 'ea23b' 'ea24b' 'ea25b' 'ea26b' 'ea27b' 'ea28b' 'ea29b' 'ea30b' ...
        'ea31b' 'ea32b' 'ea33b' 'ea34b' 'ea35b' 'ea36b' 'ea37b' 'ea38b'}; 

do_delete = 0;
cWin      = [-0.05 0.2]; % clean around blink
       
for s = 1 : length(subj)
	for bi = 1 : 4
		outname = ['ps' fs subj{s} fs subj{s} num2str(bi) '.mat'];
		if exist(outname,'file') && do_delete
			delete(outname);
		end
		if ~exist(outname,'file')
			fid = fopen(outname,'w');
			fclose(fid);
			disp(['Starting: ' outname])

			% load behavior
			log  = load(['logs' fs subj{s} num2str(bi) '.mat']);
			info = log.RT';
            no = log.no;
            
			
			resp = [];
			for ii = 1 : length(log.resp)
				if isempty(log.resp{ii})
					resp(ii) = NaN;
				elseif ismember(log.resp{ii},{'z' 'x' 'c'})
					resp(ii) = 1;
				elseif ismember(log.resp{ii},{'n' 'm' ','})
					resp(ii) = 0;
				end
			end
			info = [info resp' no];
			
			
			fname = ['rawdir' fs subj{s} num2str(bi) '_sg_res.asc'];
			[samples x y pupil resx resy tsamp tevent Sf] = read_eyelink_data_CN(fname);
%            
			ixnan = eye_clean_blinks(x,y,Sf,cWin,0);
            
			[x(ixnan) y(ixnan) pupil(ixnan) resx(ixnan) resy(ixnan)] = deal(NaN);
			% 1:samp trial onset, 2:samp fix onset, 3:con 1-noise, 0-clean, 4:trial number, 5:related 1-relat, 0-unrel, 6:duration, 7: RT, 8: resp 1-rel, 0-unrel  
			info = [tsamp(1:end,:)' tevent(1:end,:) info];
                
                
			data         = [];
			data.Sf      = Sf;
			data.time    = samples/1000;
			data.x       = x;
			data.y       = y;
			data.resx    = resx;
			data.resy    = resy;
			data.pupil   = pupil;
			data.samples = samples;
			data.events  = info;
			data.cWin    = cWin;
			save(outname,'-struct','data')
		end
	end
end
