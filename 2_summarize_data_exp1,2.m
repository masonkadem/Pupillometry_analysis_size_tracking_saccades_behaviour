%Step 3: create matrix of data from Do_read_data script (.mat file)
clear all;
format compact
fs         = filesep;

dir_data = 'PATH';
addpath(genpath([dir_data 'Exp1']))
addpath(genpath([dir_data 'Exp2']))
addpath(genpath([dir_data fs 'matlab']));
addpath(genpath([dir_data fs 'matlab' fs 'analysis']));
cd(dir_data)

% Exp1
subj = {        'ea02a'         'ea04a' 'ea05a' 'ea06a' 'ea07a' 'ea08a' 'ea09a' 'ea10a'  ...
        'ea11a' 'ea12a' 'ea13a' 'ea14a' 'ea15a' 'ea16a' 'ea17a' 'ea18a' 'ea19a' 'ea20a' ...
        'ea21a' 'ea22a' 'ea23a' 'ea24a' 'ea25a' 'ea26a'         'ea28a' 'ea29a' 'ea30a' ...
        'ea31a' 'ea32a' 'ea33a' 'ea34a' 'ea35a' 'ea36a' 'ea37a' 'ea38a'}; % excluded 27 - computer issue

% Exp2

subj = {'ea01b' 'ea02b' 'ea03b' 'ea04b' 'ea05b' 'ea06b' 'ea07b' 'ea08b' 'ea09b' 'ea10b' ...
        'ea11b' 'ea12b' 'ea13b' 'ea14b' 'ea15b' 'ea16b' 'ea17b' 'ea18b' 'ea19b' 'ea20b' ...
        'ea21b' 'ea22b' 'ea23b' 'ea24b' 'ea25b' 'ea26b' 'ea27b' 'ea28b' 'ea29b' 'ea30b' ...
        'ea31b' 'ea32b' 'ea33b' 'ea34b' 'ea35b' 'ea36b' 'ea37b' 'ea38b'}; 

tw         = [-3 7];
tshift     = 3; % time between fixation and sound onset    
do_delete = 1;

for s = 1 : length(subj)
    outname = ['ps' fs subj{s} fs subj{s} 'raw.mat'];
    
    if exist(outname,'file') & do_delete
			delete(outname);
    end
    
    disp(['Starting: ' outname])
    
    
    nn = 1;
    [epochsP events] = deal([]);
    
    for b = 1 : 4
        
        %load data
        data = load(['ps' fs subj{s} fs subj{s} num2str(b) '.mat']);
       
        for i = 1 : size(data.events(:,1),1)
            [~,ix] = min(abs(data.samples-data.events(i,1)));
                        
            tmp      = tw*data.Sf + tshift*data.Sf;
            tsamples = tmp(1) : 1 : tmp(2);
                        
            epochsP(nn,:) = data.pupil(ix+tsamples);
            epochsX(nn,:) = data.x(ix+tsamples);
            epochsY(nn,:) = data.y(ix+tsamples);
            
            % time vector for epoch
            t = tsamples/data.Sf-tshift;
            
            events(nn,:) = data.events(i,:);
            nn = nn + 1;
        end
    end
    ixtmp = isnan(events(:,9));
    events(ixtmp,8) = NaN;
    
    results = [];
    results.Sf = data.Sf;
    results.t = t;
    results.epochsP = epochsP;
    results.epochsX = epochsX;
    results.epochsY = epochsY;
    results.events = events;
    save(outname,'-struct','results')
end



