% Do_analysis_and_plotting_CN

clear all;
close all;
format compact
fs = filesep;

dir_data = 'PATH';
addpath(genpath([dir_data 'Exp1']))
addpath(genpath([dir_data 'Exp2']))
addpath(genpath([dir_data fs 'matlab']));
addpath(genpath([dir_data fs 'matlab' fs 'analysis']));
cd(dir_data)

% Exp1
% subj = {        'ea02a'         'ea04a' 'ea05a' 'ea06a' 'ea07a' 'ea08a' 'ea09a' 'ea10a'  ...
%         'ea11a' 'ea12a' 'ea13a' 'ea14a' 'ea15a' 'ea16a' 'ea17a' 'ea18a' 'ea19a' 'ea20a' ...
%         'ea21a' 'ea22a' 'ea23a' 'ea24a' 'ea25a' 'ea26a'         'ea28a' 'ea29a' 'ea30a' ...
%         'ea31a' 'ea32a' 'ea33a' 'ea34a' 'ea35a' 'ea36a' 'ea37a' 'ea38a'};

% Exp2

subj = {'ea01b' 'ea02b' 'ea03b' 'ea04b' 'ea05b' 'ea06b' 'ea07b' 'ea08b' 'ea09b' 'ea10b' ...
        'ea11b' 'ea12b' 'ea13b' 'ea14b' 'ea15b' 'ea16b' 'ea17b' 'ea18b' 'ea19b' 'ea20b' ...
        'ea21b' 'ea22b' 'ea23b' 'ea24b' 'ea25b' 'ea26b' 'ea27b' 'ea28b' 'ea29b' 'ea30b' ...
        'ea31b' 'ea32b' 'ea33b' 'ea34b' 'ea35b' 'ea36b' 'ea37b' 'ea38b'}; 


preject = 0.4;     % if 40% of trial is missing, remove
btw     = [-.5 0]; % baseline time window
tw      = [0.5 1]; % Analysis time window: first entry refers to abs time; second entry relative to sentence offset
epoch   = [-0.53 6]; % the 0.03 is because of the small edge artifact in the MS data
gwidth  = 0.02;    % ms analysis
 
% apply low pass filter
[B] = fir_design(500, 200, 'low', @kaiser, 4, 10, [0 20], [-120 10], 0.01);

[SubjDp SubjDm] = deal({});
[Dt, Dtm, Dtp, Mt, Mtm, RejData, histD] = deal([]);
for s = 1 : length(subj)
	% load data
    data1          = load(['Exp2' fs 'ps' fs subj{s} fs subj{s} 'raw' '.mat']);
    ibt            = data1.t >= epoch(1) & data1.t <= epoch(2); 
    
    % get data sorted
    Sf             = data1.Sf; % sampling frequency
    ttime          = data1.t(ibt); % time vector
    dur            = data1.events(:,7)/1000;
    condamb        = data1.events(:,2);
    condnoise      = data1.events(:,5);
    P              = data1.epochsP(:,ibt);
    X              = data1.epochsX(:,ibt);
    Y              = data1.epochsY(:,ibt);
    rt             = data1.events(:,8);  % reaction time
    isC            = double(data1.events(:,4) == data1.events(:,9)); %correct
        
    rB = [btw(1) 5];
    nB = 30;
    width = diff(rB)/nB;
    [histD(:,s) xb xbins] = linear_binning_1D(dur,rB,width,nB,@length);
    
    % calculate overal rejected data points
    RejData(s)  = nnz(isnan(P))/numel(P);
    
    % pre-process
    [PMAT XMAT YMAT] = deal(NaN(size(P)));
    for ii = 1: size(P,1)
        % remove outliers        
        out        = isoutlier(P(ii,:));
        nPfo       = P(ii,:);
        nPfo(out)  = NaN;
        
        % linearly interpolate missing values
        nPfo = fillmissing(nPfo,'linear', 'EndValues', 'nearest');
        XMAT(ii,:) = fillmissing(X(ii,:),'linear', 'EndValues', 'nearest');
        YMAT(ii,:) = fillmissing(Y(ii,:),'linear', 'EndValues', 'nearest');
        
        % filter data
        PMAT(ii,:) = filtfilt(B,1,nPfo') ; % filtering is done in first dimension
    end
    
    % get rejected trials
    [pnan pnanfull] = deal(NaN(size(dur)));
    for ii = 1 : length(dur)
		ixt = ttime >= tw(1) & ttime <= dur(ii)+tw(2);
        pnan(ii) = sum(isnan(P(ii,ixt)))/nnz(ixt);  %proportion of nan within time window
        pnanfull(ii) = sum(isnan(P(ii,:)))/size(PMAT,2);
    end

    % check whether there are some remaining NANs
    if ~isequal(sum(isnan(PMAT),2)/size(PMAT,2), pnanfull==1)
        error('Not cool')
    end
    
    % get x and y eye movement data in correct format
    M = {};
    for ii = 1 : size(XMAT,1)
        M{ii} = [XMAT(ii,:); YMAT(ii,:)];
    end
    
    % get microsaccades
    yG = [];
    MS = eye_microsaccade_detection(M,Sf,0.006,[1 1 0 -1 -1].*(Sf/6),15,1);
    for ii = 1 : length(MS)
        [y x1]    = impulses(diff(epoch),MS{ii}(:,3)/Sf,ones([length(MS{ii}(:,3)) 1]),Sf);
        yG(ii,:) = impulse_conv(y,gwidth,Sf);
    end
    x1 = x1 + epoch(1);
    
	% do baseline correction
    ixtb      = ttime >= btw(1) & ttime <= btw(2);
    D         = bsxfun(@minus,PMAT,mean(PMAT(:,ixtb),2));
    
    %rejected trials
    ixkeep    = pnan <= preject;
    ixreject  = pnan > preject;
    reject(s) = sum(ixreject)/112 ;         % proportion of trials rejected for each participant
    trialsrejected(s)   = sum(ixreject);
    totaltrialsrejected = sum(trialsrejected);
    validtrials(s)  = mean(pnan(ixkeep))  ; %mean proportion for each person
    vtm = mean(validtrials);                 % mean and 
    vts = std(validtrials);                 %STD across participants

    
    %Time course for sentences between the duration of l. NOT USED.
%     l = dur >=2 & dur <=3;
%     D23 = D(l,:);
% 
%      for ii = 0 : 1 % 0 - unamb, 1 - amb
% 		for jj = 0 : 1 % 0 - clear, 1 - noise
% 			ix = condamb(l) == ii & condnoise(l) == jj & ixkeep(l);
%             D23m(:,ii+1,jj+1,s) = mean(D23(ix,:));
%         end
%      end
    
    
    % get duration-specific mean pupil dilation
	[Dm Dp Mm] = deal(NaN(size(dur)));
    for ii = 1 : length(dur)
		ixt = ttime >= tw(1) & ttime <= dur(ii)+tw(2);
       
		Dm(ii) = mean(D(ii,ixt),2);
        Dp(ii) = max(D(ii,ixt),[],2);
        Mm(ii) = mean(yG(ii,ixt),2);
        lat    = stepinfo(D(ii,ixt), ttime(ixt));
        Pt(ii) = lat.PeakTime;
    end
    
    SubjDm{s} = Dm(ixkeep);
    SubjDp{s} = Dp(ixkeep);
    
    % get nan
    ixnan = ~isnan(data1.events(:,9)); 
    irtnan = ~isnan(data1.events(:,8)); 
    
	% average across trials: [LAC LAN; HAC HAN]
	for ii = 0 : 1 % 0 - unamb, 1 - amb
		for jj = 0 : 1 % 0 - clear, 1 - noise
			ix = condamb == ii & condnoise == jj & ixkeep;
            iy = condamb == ii & condnoise == jj & ixnan;
            ir = condamb == ii & condnoise == jj & irtnan;
            
			Dt(:,ii+1,jj+1,s) = mean(D(ix,:));
			Dtm(ii+1,jj+1,s)  = mean(Dm(ix));
            Dtp(ii+1,jj+1,s)  = mean(Dp(ix));
            
            Mt(:,ii+1,jj+1,s) = mean(yG(ix,:));
            Mtm(ii+1,jj+1,s)  = mean(Mm(ix));
            
            Ptm(ii+1,jj+1,s)  = mean(Pt(ix));   % latency
            
            rtm(ii+1,jj+1,s)  = mean(rt(ir));   % reaction time
            isCm(ii+1,jj+1,s) = mean(isC(iy));  % correct responses
            
		end
	end
end

%% Results

%Behaviour
fid = fopen('Exp1_beh_may7.txt','w');
fprintf(fid,'subjID Exp LACb HACb LANb HANb\n');
for s = 1 : length(subj)
    fprintf(fid,'%s %0d %f %f %f %f\n',subj{s},strcmp(subj{s}(5),'b')+1,isCm(1,1,s),isCm(2,1,s),isCm(1,2,s),isCm(2,2,s))
end
fclose(fid);

%pupil
fid = fopen('Exp1_pupil_may7.txt','w');
fprintf(fid,'subjID Exp LACm HACm LANm HANm LACp HACp LANp HANp LACl HACl LANl HANl\n'); 
for s = 1 : length(subj)
    fprintf(fid,'%s %0d %f %f %f %f %f %f %f %f %f %f %f %f\n',subj{s},strcmp(subj{s}(5),'b')+1,Dtm(1,1,s),Dtm(2,1,s),Dtm(1,2,s),Dtm(2,2,s),Dtp(1,1,s),Dtp(2,1,s),Dtp(1,2,s),Dtp(2,2,s), Ptm(1,1,s),Ptm(2,1,s),Ptm(1,2,s),Ptm(2,2,s));
end
fclose(fid);

%reject
fid = fopen('preject_Exp1.txt','w');
fprintf(fid,'subjdID Exp reject\n');
for s = 1: length(subj)
    fprintf(fid, '%s %0d %f\n',subj{s},strcmp(subj{s}(5),'b')+1,reject(1,s))
end
fclose(fid);


%% Stats

% do rmANOVA, post hoc and eta


% Peak latency
results = rmANOVA2(Ptm)
[~,p,~,tstats] = ttest(Ptm(1,1,:), Ptm(1,2,:)) % uc un
[~,p,~,tstats] = ttest(Ptm(2,1,:), Ptm(2,2,:)) % ac an


[~,p,~,tstats] = ttest(Ptm(1,1,:), Ptm(2,1,:)) % uc vs ac
[~,p,~,tstats] = ttest(Ptm(1,2,:), Ptm(2,2,:)) % un an

% mean
results = rmANOVA2(Dtm)
[~,p,~,tstats] = ttest(Dtm(1,1,:), Dtm(1,2,:)) % uc un
[~,p,~,tstats] = ttest(Dtm(2,1,:), Dtm(2,2,:)) % ac an


[~,p,~,tstats] = ttest(Dtm(1,1,:), Dtm(2,1,:)) % uc vs ac
[~,p,~,tstats] = ttest(Dtm(1,2,:), Dtm(2,2,:)) % un an

% peak
results = rmANOVA2(Dtp)
[~,p,~,tstats] = ttest(Dtp(1,1,:), Dtp(1,2,:)) % uc un
[~,p,~,tstats] = ttest(Dtp(2,1,:), Dtp(2,2,:)) % ac an


[~,p,~,tstats] = ttest(Dtp(1,1,:), Dtp(2,1,:)) % uc vs ac
[~,p,~,tstats] = ttest(Dtp(1,2,:), Dtp(2,2,:)) % un an
results = rmANOVA2(Mtm)

% beh
results = rmANOVA2(isCm)
[~,p,~,tstats] = ttest(isCm(1,1,:), isCm(1,2,:)) % uc un
[~,p,~,tstats] = ttest(isCm(2,1,:), isCm(2,2,:))


% eta
Noise_eta = results.SS.B/ (results.SS.B + results.SS.BxError)
Amb_eta   = results.SS.A/ (results.SS.A + results.SS.AxError)
Int_eta   = results.SS.AxB/ (results.SS.AxB + results.SS.AxBxError)

[g_eta2 c_eta2 p_eta2] = rmANOVA_eta2(SSeffect, SSerror, SSbetw)


% AU TO MM  (pupil*5)/5570.29)


% CI
% ts = tinv([0.975],length(subj)-1); 



%% % ready plot parameters

CONS = {'LAC' 'LAN'; 'HAC' 'HAN'};

% set style parameters
linstyle = {'-' ':'};
ccol     = [20 111 204; 179 9 54]/255;

% plot behaviour

% plot bar graph
x = 1:4; nn = 1; z = 1.5;
figure, set(gcf,'Color',[1 1 1]), hold on, tmpCON = {};
set(gca,'FontSize',12)
for jj = 1 : 2
    for ii = 1 : 2
		z0 = zeros(length(4), 1);
        z0(:) = z;
        v = mean(isCm(ii,jj,:));
        scatter(z0,v)
        
        errorbar(z,mean(isCm(ii,jj,:),3),sem(isCm(ii,jj,:),3),'.','Color',ccol(jj,:),'LineWidth',2)
		tmpCON{nn} = CONS{ii,jj};
        z  = z + 1;
        nn = nn + 1;
	end
end
xlim([0.8 5.2])
ylim([0.5 1])
set(gca,'XTick',x+0.5,'XTickLabel',tmpCON)
ylabel('% correct')


% plot data time course
nn = 1;
figure, set(gcf,'Color',[1 1 1]), hold on, pp = []; tmpCON = {};
set(gca,'FontSize',12)
% imagesc(xb,[-50 -45],mean(histD,2)') % default [-50 -40]
% colormap(flipud(gray(201)))
for jj = 1 : 2
    for ii = 1 : 2
		pp(nn) = plot(ttime,mean(D23m(:,ii,jj,:),4),linstyle{ii},'Color',ccol(jj,:),'LineWidth',2);
		tmpCON{nn} = CONS{ii,jj};
		nn = nn + 1;
	end
end
xlim([btw(1) 4])
xlabel('Time (s)')
ylabel('Pupil dilation (a.u.)')
legend(pp,tmpCON,'Location','NorthWest')
legend('boxoff')
%colorbar


% plot bar graph
x = 1:4; nn = 1;
figure, set(gcf,'Color',[1 1 1]), hold on, tmpCON = {};
set(gca,'FontSize',12)
for jj = 1 : 2
    for ii = 1 : 2
		rectangle('Position',[x(nn) 0 1 mean(Dtm(ii,jj,:),3)],'FaceColor',[0.8 0.8 0.8],'EdgeColor',ccol(jj,:),'LineStyle',linstyle{ii},'LineWidth',2)
		errorbar(x(nn)+0.5,mean(Dtm(ii,jj,:),3),sem(Dtm(ii,jj,:),3),'.','Color',ccol(jj,:),'LineWidth',2)
		tmpCON{nn} = CONS{ii,jj};
		nn = nn + 1;
	end
end
xlim([0.8 5.2])
set(gca,'XTick',x+0.5,'XTickLabel',tmpCON)
ylabel('Pupil dilation (a.u.)')


% plot bar peak graph
x = 1:4; nn = 1;
figure, set(gcf,'Color',[1 1 1]), hold on, tmpCON = {};
set(gca,'FontSize',12)
for jj = 1 : 2
    for ii = 1 : 2
		rectangle('Position',[x(nn) 0 1 mean(Dtp(ii,jj,:),3)],'FaceColor',[0.8 0.8 0.8],'EdgeColor',ccol(jj,:),'LineStyle',linstyle{ii},'LineWidth',2)
		errorbar(x(nn)+0.5,mean(Dtp(ii,jj,:),3),sem(Dtp(ii,jj,:),3),'.','Color',ccol(jj,:),'LineWidth',2)
		tmpCON{nn} = CONS{ii,jj};
		nn = nn + 1;
	end
end
xlim([0.8 5.2])
set(gca,'XTick',x+0.5,'XTickLabel',tmpCON)
ylabel('Peak pupil dilation (a.u.)')


% plot scatter
% Ambiguity
x = 1:4; nn = 1;
figure, set(gcf,'Color',[1 1 1]), hold on, tmpCON = {};
set(gca,'FontSize',12)
subplot(1,2,1), hold on
plot([-200 400],[-200 400],':k')
scatter(mean(Dtm(1,:,:),2), mean(Dtm(2,:,:),2),80,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0]);
ylabel('High ambiguity ')
xlabel('Low ambiguity ')
daspect([1 1 1])


% Noise
subplot(1,2,2), hold on
plot([-200 400],[-200 400],':k')
scatter(mean(Dtm(:,1,:),1), mean(Dtm(:,2,:),1),80,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0]);
ylabel('Noise')
xlabel('Clear')
daspect([1 1 1])


% plot micro saccade time course
nn = 1;
figure, set(gcf,'Color',[1 1 1]), hold on, pp = []; tmpCON = {};
set(gca,'FontSize',12)
for jj = 1 : 2
    for ii = 1 : 2
		pp(nn) = plot(x1,mean(Mt(:,ii,jj,:),4),linstyle{ii},'Color',ccol(jj,:),'LineWidth',2);
		tmpCON{nn} = CONS{ii,jj};
		nn = nn + 1;
	end
end
xlim([btw(1) 4])
xlabel('Time (s)')
ylabel('Frequency')
legend(pp,tmpCON,'Location','NorthWest')
legend('boxoff')


% plot bar graph
x = 1:4; nn = 1;
figure, set(gcf,'Color',[1 1 1]), hold on, tmpCON = {};
set(gca,'FontSize',12)
for jj = 1 : 2
    for ii = 1 : 2
		rectangle('Position',[x(nn) 0 1 mean(Mtm(ii,jj,:),3)],'FaceColor',[0.8 0.8 0.8],'EdgeColor',ccol(jj,:),'LineStyle',linstyle{ii},'LineWidth',2)
		errorbar(x(nn)+0.5,mean(Mtm(ii,jj,:),3),sem(Mtm(ii,jj,:),3),'.','Color',ccol(jj,:),'LineWidth',2)
		tmpCON{nn} = CONS{ii,jj};
		nn = nn + 1;
	end
end
xlim([0.8 5.2])
set(gca,'XTick',x+0.5,'XTickLabel',tmpCON)
ylabel('Frequency')


% scatter plot, mean vs. peak, 3 participants



% fid = fopen('filename1.txt','w');
% fprintf(fid,'subjID Exp UCm ACm UNm ANm UCp ACp UNp ANp\n');
% for s = 1 : length(subj)
%     fprintf(fid,'%s %0d %f %f %f %f %f %f %f %f\n',subj{s},strcmp(subj{s}(5),'b')+1,Dtm(1,1,s),Dtm(2,1,s),Dtm(1,2,s),Dtm(2,2,s),Dtp(1,1,s),Dtp(2,1,s),Dtp(1,2,s),Dtp(2,2,s));
% end
% fclose(fid);


