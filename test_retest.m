clear all; close all; clc

WhichParc = 'Gordon' % 'Gordon' 'Power'
switch WhichParc
	case 'Gordon'
		Parc = 1;
	case 'Power'
		Parc = 2;
end

RunGroupDiff = false
WhichStat = 't'

% NOTE: only run one or the other. Both CANNOT be set to true
runSR = false
runScrub = false

% ------------------------------------------------------------------------------
% This script performs test-retest analysis on the NYU_2 dataset
% 
% Linden Parkes, Brain & Mental Health Laboratory, 2016
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Project
% ------------------------------------------------------------------------------
sublist = '~/Dropbox/Work/Scripts/Sublists/NYU_2.txt';
projdir = '~/Dropbox/Work/ResProjects/2016/rfMRI/NYU_2/';
datadir = [projdir,'data/'];
% Baseline data directory string
% Note, we use the baseline data to calculate motion
preprostr = '/session_1/rest_1/prepro/';

% ------------------------------------------------------------------------------
% Subject list
% ------------------------------------------------------------------------------
fileID = fopen(sublist);
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};

% compute numsubs
numSubs = length(ParticipantIDs);

% ------------------------------------------------------------------------------
% Jenkinson's mean FD
% ------------------------------------------------------------------------------
fprintf(1, 'Loading Jenkinson''s mean FD metric\n');
[exclude,~,fdJenk,fdJenk_m] = RunExclude(datadir,ParticipantIDs,preprostr);
fprintf(1, 'done\n');

% Sort movement ascending
[fdJenk_m_srt,idx] = sort(fdJenk_m);

% rearrange daris IDs
ParticipantIDs = ParticipantIDs(idx);

% create grouping variable based on motion
% 1 = low, 2 = medium, 3 = high
x = round(numSubs/3);
Group = [zeros(x,1)+1; zeros(x,1)+2; zeros(x,1)+3];
if length(Group) > numSubs
	Group(numSubs+1:end) = [];
elseif length(Group) < numSubs
	Group(numSubs) = 3;
	Group(end-10:end) = 3;
end

% ------------------------------------------------------------------------------
% Load ROI coordinates
% ------------------------------------------------------------------------------
switch WhichParc
	case 'Gordon'
		% ROI_Coords = dlmread('/gpfs/M2Home/projects/Monash076/Linden/ROIs/Gordon/Gordon_Centroids.txt')
		ROI_Coords = dlmread('~/Dropbox/Work/ROIs/Gordon/Gordon_Centroids.txt');
	case 'Power'
		ROI_Coords = dlmread('~/Dropbox/Work/ROIs/Power/Power2011_xyz_711.txt');
end

% Calculate pairwise euclidean distance
ROIDist = GetROIDist(ROI_Coords);

% Flatten distance matrix
ROIDistVec = LP_FlatMat(ROIDist);

% Calculate number of ROIs
numROIs = size(ROIDist,1);

% Calculate number of edges
numConnections = numROIs * (numROIs - 1) / 2;

% ------------------------------------------------------------------------------
% Noise corrections
% ------------------------------------------------------------------------------
	noiseOptions = {'6P',...
					'6P+2P',...
					'6P+2P+GSR',...
					'24P',...
					'24P+8P',...
					'24P+8P+4GSR',...
					'24P+aCC',...
					'24P+aCC50',...
					'24P+aCC+4GSR',...
					'24P+aCC50+4GSR',...
					'12P+aCC',...
					'12P+aCC50',...
					'sICA-AROMA+2P',...
					'sICA-AROMA+2P+GSR',...
					'sICA-AROMA+8P',...
					'sICA-AROMA+8P+4GSR'};
	noiseOptionsNames = {'6HMP',...
						'6HMP+2Phys',...
						'6HMP+2Phys+GSR',...
						'24HMP',...
						'24HMP+8Phys',...
						'24HMP+8Phys+4GSR',...
						'24HMP+aCompCor',...
						'24HMP+aCompCor50',...
						'24HMP+aCompCor+4GSR',...
						'24HMP+aCompCor50+4GSR',...
						'12HMP+aCompCor',...
						'12HMP+aCompCor50',...
						'ICA-AROMA+2Phys',...
						'ICA-AROMA+2Phys+GSR',...
						'ICA-AROMA+8Phys',...
						'ICA-AROMA+8Phys+4GSR'};

numPrePro = length(noiseOptions);

cfgFile = 'cfg.mat';

% ------------------------------------------------------------------------------
% Variables
% ------------------------------------------------------------------------------
allStats = cell(numPrePro,2);
allStats = struct('ICC_wrt',cell(1,numPrePro),...
					'ICC_brt',cell(1,numPrePro),...
					'ICC_wrt_mean',[],...
					'ICC_wrt_std',[],...
					'ICC_brt_mean',[],...
					'ICC_brt_std',[],...
					'tstat',[],...
					'p',[],...
					'p_corrected',[],...
					'tDOF_mean',[],...
					'tDOF_std',[]);

% ------------------------------------------------------------------------------
% Figures
% ------------------------------------------------------------------------------
% Figure 1
figure('color','w', 'units', 'centimeters', 'pos', [0 0 50 10], 'name',['']); box('on'); hold on;

% ------------------------------------------------------------------------------
% Script
% ------------------------------------------------------------------------------
for i = 1:numPrePro
    removeNoise = noiseOptions{i};
	fprintf(1, '\nProcessing data: %s\n',removeNoise);

	% ------------------------------------------------------------------------------
	% 1) Baseline data
	% ------------------------------------------------------------------------------
	    cfg_bl = [];
	    % preprostr
	    preprostr = '/session_1/rest_1/prepro/';

		% ------------------------------------------------------------------------------
		% Load in time series data
		% ------------------------------------------------------------------------------
		for j = 1:numSubs
		    tsdir = [datadir,ParticipantIDs{j},preprostr,removeNoise,'/'];
		    
		    clear temp
		    temp = load([tsdir,cfgFile]);
		    
		    cfg_bl = [cfg_bl temp.cfg];
		end

	    % ------------------------------------------------------------------------------
	    % Compute correlations
	    % ------------------------------------------------------------------------------
		numROIs = size(cfg_bl(1).roiTS{Parc},2);

		FC_bl = zeros(numROIs,numROIs,numSubs);
		FC_bl_vec = zeros(numSubs,numConnections);

		for j = 1:numSubs
			FC_bl(:,:,j) = corr(cfg_bl(j).roiTS{Parc});
			% Perform fisher z transform
			FC_bl(:,:,j) = fisherz(FC_bl(:,:,j));
			% Flatten
			FC_bl_vec(j,:) = LP_FlatMat(FC_bl(:,:,j));
		end

	% ------------------------------------------------------------------------------
	% Get tDOF
	% ------------------------------------------------------------------------------
	% we use baseline data to compute tDOF

	fprintf(1, 'Computing tDOF: %s\n',removeNoise);
	tDOFtemp = zeros(numSubs,1);
	for j = 1:numSubs

		% get tDOF
		% First, find size of second dimension of noiseTS
		if runSR
			tDOFtemp(j) = size(cfg_bl(j).noiseTS_spikereg,2);
		else
			tDOFtemp(j) = size(cfg_bl(j).noiseTS,2);
		end

		if runScrub
			tDOFtemp(j) = tDOFtemp(j) + sum(ScrubMask{j});
		end

		% Then, if ICA-AROMA pipeline, find number of ICs and add to tDOF
		if ~isempty(strfind(removeNoise,'ICA-AROMA'))
			if runSR
				x = dlmread([datadir,ParticipantIDs{j},preprostr,removeNoise,'+SpikeReg/classified_motion_ICs.txt']);
			else
				x = dlmread([datadir,ParticipantIDs{j},preprostr,removeNoise,'/classified_motion_ICs.txt']);
			end
			tDOFtemp(j) = tDOFtemp(j) + length(x);
		end
	end

	% ------------------------------------------------------------------------------
	% Calculate mean temporal degrees of freedom lost
	% ------------------------------------------------------------------------------
	% tDOF will be the same for most pipelines, but some have variable regressor amounts
	% So we take mean over subjects
	allStats(i).tDOF_mean = mean(tDOFtemp);
	allStats(i).tDOF_std = std(tDOFtemp);

	% ------------------------------------------------------------------------------
	% 2) Within session retest data
	% ------------------------------------------------------------------------------
	    cfg_wrt = [];
	    % preprostr
	    preprostr = '/session_1/rest_2/prepro/';

		% ------------------------------------------------------------------------------
		% Load in time series data
		% ------------------------------------------------------------------------------
		for j = 1:numSubs
		    tsdir = [datadir,ParticipantIDs{j},preprostr,removeNoise,'/'];
		    
		    clear temp
		    temp = load([tsdir,cfgFile]);
		    
		    cfg_wrt = [cfg_wrt temp.cfg];
		end

	    % ------------------------------------------------------------------------------
	    % Compute correlations
	    % ------------------------------------------------------------------------------
		numROIs = size(cfg_wrt(1).roiTS{Parc},2);

		FC_wrt = zeros(numROIs,numROIs,numSubs);
		FC_wrt_vec = zeros(numSubs,numConnections);

		for j = 1:numSubs
			FC_wrt(:,:,j) = corr(cfg_wrt(j).roiTS{Parc});
			% Perform fisher z transform
			FC_wrt(:,:,j) = fisherz(FC_wrt(:,:,j));
			% Flatten
			FC_wrt_vec(j,:) = LP_FlatMat(FC_wrt(:,:,j));
		end

	% ------------------------------------------------------------------------------
	% 3) Between session retest data
	% ------------------------------------------------------------------------------
	    cfg_brt = [];
	    % preprostr
	    preprostr = '/session_2/rest_1/prepro/';

		% ------------------------------------------------------------------------------
		% Load in time series data
		% ------------------------------------------------------------------------------
		for j = 1:numSubs
		    tsdir = [datadir,ParticipantIDs{j},preprostr,removeNoise,'/'];
		    
		    clear temp
		    temp = load([tsdir,cfgFile]);
		    
		    cfg_brt = [cfg_brt temp.cfg];
		end

	    % ------------------------------------------------------------------------------
	    % Compute correlations
	    % ------------------------------------------------------------------------------
		numROIs = size(cfg_brt(1).roiTS{Parc},2);

		FC_brt = zeros(numROIs,numROIs,numSubs);
		FC_brt_vec = zeros(numSubs,numConnections);

		for j = 1:numSubs
			FC_brt(:,:,j) = corr(cfg_brt(j).roiTS{Parc});
			% Perform fisher z transform
			FC_brt(:,:,j) = fisherz(FC_brt(:,:,j));
			% Flatten
			FC_brt_vec(j,:) = LP_FlatMat(FC_brt(:,:,j));
		end

	% ------------------------------------------------------------------------------
	% Calculate ICC
	% ------------------------------------------------------------------------------
		fprintf(1, 'Performing ICC for: %s',removeNoise);
		% ICC_wrt = arrayfun(@(x)GetICC([FC_bl_vec(:,x),FC_wrt_vec(:,x)]),1:numConnections);

		% initialise
		allStats(i).ICC_wrt = zeros(1,numConnections);
		allStats(i).ICC_brt = zeros(1,numConnections);

		% For each functional connection
		for j = 1:numConnections
			% 1) within session T-RT
			x = [FC_bl_vec(:,j) FC_wrt_vec(:,j)];
			allStats(i).ICC_wrt(j) = GetICC(x);

			% 2) between session T-RT
			y = [FC_bl_vec(:,j) FC_brt_vec(:,j)];
			allStats(i).ICC_brt(j) = GetICC(y);
		end

		% mean
		allStats(i).ICC_wrt_mean = nanmean(allStats(i).ICC_wrt);
		allStats(i).ICC_wrt_std = nanstd(allStats(i).ICC_wrt);

		% std
		allStats(i).ICC_brt_mean = nanmean(allStats(i).ICC_brt);
		allStats(i).ICC_brt_std = nanstd(allStats(i).ICC_brt);

		fprintf(1, '...done \n');

	% ------------------------------------------------------------------------------
	% Plot over distance (Figure 1)
	% ------------------------------------------------------------------------------
		fprintf(1, 'Drawing distance plot \n');
	
		figure(1)
		subplot(2,round(numPrePro/2),i)

		numThresholds = 10;
		xData = 1:numThresholds - 1;

		[yMeans,yStds] = LP_BinData(ROIDistVec',allStats(i).ICC_wrt,numThresholds);
		boundedline(xData,yMeans,yStds,'-bo','alpha')

		[yMeans,yStds] = LP_BinData(ROIDistVec',allStats(i).ICC_brt,numThresholds);
		boundedline(xData,yMeans,yStds,'-ro','alpha')

		ylim([0 1])
		title(removeNoise,'Interpreter', 'none')


	% ------------------------------------------------------------------------------
	% Compute within session and between session ICC differences
	% ------------------------------------------------------------------------------
		switch WhichStat
			case 't'
				[~,p,~,stats] = ttest(allStats(i).ICC_wrt,allStats(i).ICC_brt);

				allStats(i).tstat = stats.tstat;
				allStats(i).p = p;
		end

	% ------------------------------------------------------------------------------
	% Calculate high/low motion differences
	% Note, we use the baseline data to calculate motion
	% ------------------------------------------------------------------------------
	if RunGroupDiff
		fprintf(1, 'Calculating motion group differences \n');
		numGroups = numel(unique(Group));

		% Calculate ICC on low and high motion groups separately

		% initialise
		allStats(i).ICC_wrt_mot = cell(1,numGroups);
		allStats(i).ICC_brt_mot = cell(1,numGroups);

		% For each motion group
		for g = [1,numGroups]
			% initialise
			allStats(i).ICC_wrt_mot{g} = zeros(1,numConnections);
			allStats(i).ICC_brt_mot{g} = zeros(1,numConnections);

			% For each functional connection
			for j = 1:numConnections
				% 1) within session T-RT
				x = [FC_bl_vec(Group == g,j) FC_wrt_vec(Group == g,j)];
				allStats(i).ICC_wrt_mot{g}(j) = GetICC(x);

				% 2) between session T-RT
				y = [FC_bl_vec(Group == g,j) FC_brt_vec(Group == g,j)];
				allStats(i).ICC_brt_mot{g}(j) = GetICC(y);
			end
		end		

		% calculate test stat for ICC between high/low motion
		for j = 1:2
			if j == 1
				% 1) WRT ICC
				x = allStats(i).ICC_wrt_mot{1}; % low
				y = allStats(i).ICC_wrt_mot{3}; % high
			elseif j == 2
				% 2) BRT ICC
				x = allStats(i).ICC_brt_mot{1}; % low
				y = allStats(i).ICC_brt_mot{3}; % high
			end

			switch WhichStat
				case 't'
					[~,p,~,stats] = ttest2(x,y,'Vartype','unequal');

					allStats(i).tstat_mot(j) = stats.tstat;
					allStats(i).p_mot(j) = p;
				case 'u'
					% ranksum requires inputs be a vector so we have to loop over variables (columns)
					% Loop over edges (columns)
					for j = 1:size(x,2)
						[p,~,stats] = ranksum(x(:,j),y(:,j));
						% Normalized Mann-Whitney U test (given the sample size may change across features)
						n1 = sum(Group == g1); n2 = sum(Group == g2);
						uStat = stats.ranksum - n1*(n1+1)/2;
						
						allStats{i,j}(j,1) = uStat/n1/n2;
						allStats{i,j}(j,2) = p;
					end
					% correct p values using FDR
					P_corrected = mafdr(allStats{i,j}(:,2),'BHFDR','true');
					allStats{i,j}(:,3) = P_corrected;
				end
		end
	end
end

% ------------------------------------------------------------------------------
% Correct t-stats for multiple comparisons
% ------------------------------------------------------------------------------
	% correct p values using FDR
	x = mafdr([allStats(:).p],'BHFDR','true');
	for i = 1:numPrePro
		allStats(i).p_corrected = x(i);
	end

	if RunGroupDiff
		Y = vertcat(allStats(:).p_mot);
		y = [mafdr(Y(:,1),'BHFDR','true') mafdr(Y(:,2),'BHFDR','true')];
		for i = 1:numPrePro
			allStats(i).p_mot_corrected = y(i,:);
		end
	end

% ------------------------------------------------------------------------------
% Violin plot (Figure 2)
% ------------------------------------------------------------------------------
	% Figure 2
	figure('color','w', 'units', 'centimeters', 'pos', [0 0 50 10], 'name',['']); box('on'); hold on;
	
	extraParams.theColors = [0 0.8 0; 0 0 0.8];
	extraParams.theColors = repmat(extraParams.theColors,numPrePro,1);
	extraParams.theColors = num2cell(extraParams.theColors,2);

	extraParams.customSpot = '';

	extraParams.theLabels = reshape(repmat(noiseOptions,2,1),1,numPrePro*2);

	% Build Y
	Y = reshape({allStats(:).ICC_wrt;allStats(:).ICC_brt},1,numPrePro*2);

	starPositions = [1:2:numPrePro*2]+0.5;

	JitteredParallelScatter(Y,1,1,0,extraParams)
	set(gca,'FontSize',15,'fontWeight','bold')
	title('Within and between session ICC','Interpreter', 'none')

	% Add stars if high/low ICC is significantly different
	for i = 1:numPrePro
        if allStats(i).p_corrected < 0.05
            line([starPositions(i) - 0.25,starPositions(i) + 0.25],[0.9,0.9],'Color','k'); hold on
            scatter(starPositions(i),0.8,'*','k'); hold on
        end
    end

% ------------------------------------------------------------------------------
% Sort tDOF
% ------------------------------------------------------------------------------
[~,tDOF_idx] = sort([allStats.tDOF_mean],'ascend');

% ------------------------------------------------------------------------------
% WRT ICC
% ------------------------------------------------------------------------------
	% Create data
	data = {[allStats(:).ICC_wrt_mean]'};
	data_std = {[allStats(:).ICC_wrt_std]'};

	% reorder by tDOF-loss
	data{1} = data{1}(tDOF_idx);
	data_std{1} = data_std{1}(tDOF_idx);

	% Create table
	T = table(data{1},'RowNames',noiseOptions(tDOF_idx)','VariableNames',{'WRT_ICC'})

	% Create bar chart
	clear extraParams
	extraParams.xTickLabels = noiseOptionsNames(tDOF_idx);
	extraParams.XTickLabelRot = 90;
	extraParams.xLabel = 'Pipeline';
	extraParams.yLabel = 'Within ICC';
	extraParams.xLimits = [0 numPrePro+1];
	extraParams.yLimits = [0 1];
	extraParams.Title = ['NYU dataset. ',WhichParc,' parcellation'];

	TheBarChart(data,data_std,extraParams)

% ------------------------------------------------------------------------------
% BRT ICC
% ------------------------------------------------------------------------------
	% Create data
	data = {[allStats(:).ICC_brt_mean]'};
	data_std = {[allStats(:).ICC_brt_std]'};
	
	% reorder by tDOF-loss
	data{1} = data{1}(tDOF_idx);
	data_std{1} = data_std{1}(tDOF_idx);

	% Create table
	T = table(data{1},'RowNames',noiseOptions(tDOF_idx)','VariableNames',{'BRT_ICC'})

	% Create bar chart
	clear extraParams
	extraParams.xTickLabels = noiseOptionsNames(tDOF_idx);
	extraParams.XTickLabelRot = 90;
	extraParams.xLabel = 'Pipeline';
	extraParams.yLabel = 'Between ICC';
	extraParams.xLimits = [0 numPrePro+1];
	extraParams.yLimits = [0 1];
	extraParams.Title = ['NYU dataset. ',WhichParc,' parcellation'];

	TheBarChart(data,data_std,extraParams)

% ------------------------------------------------------------------------------
% ICC tstat
% ------------------------------------------------------------------------------
	% Create data
	data = {[allStats(:).tstat]'};
	data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
	
	% reorder by tDOF-loss
	data{1} = data{1}(tDOF_idx);
	data_std{1} = data_std{1}(tDOF_idx);

	% Create table
	T = table(data{1},'RowNames',noiseOptions(tDOF_idx)','VariableNames',{'ICC_tstat'})

	% Create bar chart
	clear extraParams
	extraParams.xTickLabels = noiseOptionsNames(tDOF_idx);
	extraParams.XTickLabelRot = 90;
	extraParams.xLabel = 'Pipeline';
	extraParams.yLabel = 'ICC t-stat';
	extraParams.xLimits = [0 numPrePro+1];
	extraParams.yLimits = [0 230];
	extraParams.Title = ['NYU dataset. ',WhichParc,' parcellation'];

	TheBarChart(data,data_std,extraParams)

% ------------------------------------------------------------------------------
% Plot high/low (Figure 3)
% ------------------------------------------------------------------------------
if RunGroupDiff
	% Figure 3
	figure('color','w', 'units', 'centimeters', 'pos', [0 0 50 10], 'name',['']); box('on'); hold on;

	extraParams.theColors = [0 0.8 0; 0 0 0.8];
	extraParams.theColors = repmat(extraParams.theColors,numPrePro,1);
	extraParams.theColors = num2cell(extraParams.theColors,2);

	extraParams.customSpot = '';

	extraParams.theLabels = reshape(repmat(noiseOptions,2,1),1,numPrePro*2);
	
	% Build Y and plot
	% For each TRT (within and between)

	figureLabels = {'within TRT','between TRT'};
	
	% for each ICC
	for j = 1:2
		if j == 1
			% 1) WRT ICC
			X = vertcat(allStats(:).ICC_wrt_mot);
			x = vertcat({X{:,1}},{X{:,3}});
			Y = reshape(x,1,numPrePro*2);
		elseif j == 2
			% 2) BRT ICC
			X = vertcat(allStats(:).ICC_brt_mot);
			x = vertcat({X{:,1}},{X{:,3}});
			Y = reshape(x,1,numPrePro*2);
		end

		subplot(2,1,j)
		JitteredParallelScatter(Y,1,1,0,extraParams)
		title(figureLabels{j},'Interpreter', 'none')
		hold on

		% Add stars if high/low ICC is significantly different
		for i = 1:numPrePro
	        if allStats(i).p_mot_corrected(j) < 0.05
	            line([starPositions(i) - 0.25,starPositions(i) + 0.25],[0.9,0.9],'Color','k'); hold on
	            scatter(starPositions(i),0.8,'*','k'); hold on
	        end
	    end
	end
end	

% ------------------------------------------------------------------------------
% Create table
% ------------------------------------------------------------------------------
	varNames = {'ICC_Within_m_sd','ICC_Between_m_sd'};
	T = table([round([allStats(:).ICC_wrt_mean],2)',round([allStats(:).ICC_wrt_std],2)'],[round([allStats(:).ICC_brt_mean],2)',round([allStats(:).ICC_brt_std]',2)],'RowNames',noiseOptions','VariableNames',varNames')




