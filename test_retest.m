clear all; close all; clc
WhichProject = 'NYU_2'
% ------------------------------------------------------------------------------
% This script performs test-retest analysis on the NYU_2 dataset
% 
% Linden Parkes, Brain & Mental Health Laboratory, 2016
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Select project
% ------------------------------------------------------------------------------
switch WhichProject
	case 'NYU_2'
		sublist = '/gpfs/M2Home/projects/Monash076/Linden/Sublists/NYU_2.txt';
        datadir = '/gpfs/M2Home/projects/Monash076/Linden/NYU_2/data/';
end

% ------------------------------------------------------------------------------
% Subject list
% ------------------------------------------------------------------------------
fileID = fopen(sublist);
DarisIDs = textscan(fileID,'%s');
DarisIDs = DarisIDs{1};

% compute numsubs
numSubs = length(DarisIDs);

% ------------------------------------------------------------------------------
% Load ROI coordinates
% ------------------------------------------------------------------------------
load('/gpfs/M2Home/projects/Monash076/Linden/ROIs/Gordon/Gordon_Centroids.mat')
% load('~/Dropbox/Work/ROIs/Gordon/Gordon_Centroids.mat')

% Calculate pairwise euclidean distance
GordonDist = GetROIDist(Gordon_Coords);

% Flatten distance matrix
GordonDistVec = LP_FlatMat(GordonDist)';

% ------------------------------------------------------------------------------
% Noise corrections
% ------------------------------------------------------------------------------
% noiseOptions = {'6P','6P+2P','6P+2P+GSR','24P+8P','24P+8P+GSR','24P+aCC','24P+aCCd','24P+8P+GSR+SpikeReg','sICA-AROMA+2P','sICA-AROMA+GSR'};
noiseOptions = {'6P','6P+2P','6P+2P+GSR','24P+8P','24P+8P+GSR','24P+aCC','24P+aCCd','24P+8P+GSR+SpikeReg','ICA-AROMA+2P','ICA-AROMA+GSR'};
% noiseOptions = {'sICA-AROMA+2P','sICA-AROMA+GSR','ICA-AROMA+2P','ICA-AROMA+GSR'};
numPrePro = length(noiseOptions);
WhichParc = 1;

cfgFile = 'cfg_noSmooth.mat';

mean_wrt_ICC = [];
std_wrt_ICC = [];

mean_brt_ICC = [];
std_brt_ICC = [];

figure('color','w', 'units', 'centimeters', 'pos', [0 0 50 10], 'name',['']); box('on'); hold on;

for i = 1:numPrePro
    removeNoise = noiseOptions{i};
	fprintf(1, 'Performing ICC for: %s\n',removeNoise);

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
		    tsdir = [datadir,DarisIDs{j},preprostr,removeNoise,'/'];
		    
		    clear temp
		    temp = load([tsdir,cfgFile]);
		    
		    cfg_bl = [cfg_bl temp.cfg];
		end

	    % ------------------------------------------------------------------------------
	    % Compute correlations
	    % ------------------------------------------------------------------------------
		numROIs = size(cfg_bl(1).roiTS{WhichParc},2);

		FC_bl = zeros(numROIs,numROIs,numSubs);
		FC_bl_vec = []; 
		for j = 1:numSubs
			FC_bl(:,:,j) = corr(cfg_bl(j).roiTS{WhichParc});
			% Perform fisher z transform
			FC_bl(:,:,j) = fisherz(FC_bl(:,:,j));
			% Flatten
			FC_bl_vec(j,:) = LP_FlatMat(FC_bl(:,:,j));
		end

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
		    tsdir = [datadir,DarisIDs{j},preprostr,removeNoise,'/'];
		    
		    clear temp
		    temp = load([tsdir,cfgFile]);
		    
		    cfg_wrt = [cfg_wrt temp.cfg];
		end

	    % ------------------------------------------------------------------------------
	    % Compute correlations
	    % ------------------------------------------------------------------------------
		numROIs = size(cfg_wrt(1).roiTS{WhichParc},2);

		FC_wrt = zeros(numROIs,numROIs,numSubs);
		FC_wrt_vec = [];
		for j = 1:numSubs
			FC_wrt(:,:,j) = corr(cfg_wrt(j).roiTS{WhichParc});
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
		    tsdir = [datadir,DarisIDs{j},preprostr,removeNoise,'/'];
		    
		    clear temp
		    temp = load([tsdir,cfgFile]);
		    
		    cfg_brt = [cfg_brt temp.cfg];
		end

	    % ------------------------------------------------------------------------------
	    % Compute correlations
	    % ------------------------------------------------------------------------------
		numROIs = size(cfg_brt(1).roiTS{WhichParc},2);

		FC_brt = zeros(numROIs,numROIs,numSubs);
		FC_brt_vec = [];
		for j = 1:numSubs
			FC_brt(:,:,j) = corr(cfg_brt(j).roiTS{WhichParc});
			% Perform fisher z transform
			FC_brt(:,:,j) = fisherz(FC_brt(:,:,j));
			% Flatten
			FC_brt_vec(j,:) = LP_FlatMat(FC_brt(:,:,j));
		end

	% ------------------------------------------------------------------------------
	% Calculate ICC
	% ------------------------------------------------------------------------------
		numConnections = numROIs * (numROIs - 1) / 2;
		
		wrt_ICC = [];
		brt_ICC = [];
		% For each functional connection
		for j = 1:numConnections
			% 1) within session T-RT
			x = [FC_bl_vec(:,j) FC_wrt_vec(:,j)];
			wrt_ICC(j) = GetICC(x);

			% 2) between session T-RT
			y = [FC_bl_vec(:,j) FC_brt_vec(:,j)];
			brt_ICC(j) = GetICC(y);
		end

		% mean
		mean_wrt_ICC(i) = nanmean(wrt_ICC);
		std_wrt_ICC(i) = nanstd(wrt_ICC);

		mean_brt_ICC(i) = nanmean(brt_ICC);
		std_brt_ICC(i) = nanstd(brt_ICC);

	% ------------------------------------------------------------------------------
	% Plot
	% ------------------------------------------------------------------------------
		subplot(2,round(numPrePro/2),i)
		% histogram(wrt_ICC); hold on; histogram(brt_ICC)

		numThresholds = 10;
		xData = 1:numThresholds - 1;

		[yMeans,yStds] = LP_BinData(GordonDistVec,wrt_ICC,numThresholds);
		boundedline(xData,yMeans,yStds,'-bo','alpha')

		[yMeans,yStds] = LP_BinData(GordonDistVec,brt_ICC,numThresholds);
		boundedline(xData,yMeans,yStds,'-ro','alpha')

		ylim([0 1])
		title(removeNoise,'Interpreter', 'none')

end