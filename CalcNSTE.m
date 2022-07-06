function NSTE = CalcNSTE(Target, Source, m, tau, pred_time, nshuf)
%% CALCNSTE - normalized symbolic transfer entropy calculation
% Calculates the normalized symbolic transfer entropy for the input data 
% using the methods outlined in Lee et al., 2013, Anesthesiology. The
% current implementation only works for single trial data.
%
% Syntax: [NSTE] = CalcNSTE(Target, Source, m, tau, pred_time, nshuf)
%
% Input:
%   Target (double) - Raw data corresponding to the target channel of the
%       NSTE calculation.
%
%   Source (double) - Raw data corresponding to the source channel of the
%       NSTE calculation.
%
%   m (integer) - Integer value determining the embedding dimension of the
%       NSTE calculation, i.e. the number of data points in each vector
%       being symbolized.
%
%   tau (integer) - Integer value determining the time scale of information
%       transfer. Corresponds to the number of samples between datapoints 
%       used for symbolization. Making this larger increases the time-scale 
%       of information transfer. Default is 1 (i.e. all neighboring points 
%       used). 
%
%   pred_time (integer) - Integer value determining how far ahead to
%       predict during the transfer entropy calculation. Corresponds to the
%       'delay_rec' parameter in the quickTE function.
%
%   nshuf (integer) - The number of times the data should be shuffled to
%       estimate the bias in the data when performing the normalization
%       calculation. Increasing this number will lead to more accurate
%       estimation of bias but will increase computational times. Lee et 
%       al. used nshuf = 20.
%       
% Outputs:
%   NSTE (double) - the normalized symbolic transfer entropy from Source to
%     Target. Value is bound between [0,1].
%   
% Examples:
%   Single Trial Data
%     nste = CalcNSTE(y, x, 3, 1, 5, 20);
%
% Other m-files required: symbolicdata2states, quickTE, quickEnt
% Subfunctions: none
% MAT-files required: none
%
symbinfo = [1,m];
MethodAssign = {1,1,'Nat';1,2,'Nat'}; %Tells quickTE to not bin the data since we have already symbolized it

%transpose data if necessary
if size(Target,1) > 1; Target = Target'; end
if size(Source,1) > 1; Source = Source'; end

% Symbolic Transfer Entropy
dataRaster = [Target; Source];
symbRaster = symbolicdata2states(dataRaster,symbinfo, 'tau', tau); %First symbolize the EEG data using embedding dimension m
STE = quickTE(symbRaster(1,:)',symbRaster(2,:)','delayRec', pred_time,...
              'delayTrans',pred_time,'MethodAssign',MethodAssign);

% Symbolic Transfer Entropy of shuffled data
tmp = NaN(1,nshuf);
for n = 1:nshuf
dataRaster = [Target; Source(randperm(length(Source)))]; %Shuffle the 'Source' channel to eliminate causality
symbRaster = symbolicdata2states(dataRaster,symbinfo, 'tau', tau); %First symbolize the EEG data using embedding dimension m
tmp(n) = quickTE(symbRaster(1,:)',symbRaster(2,:)','delayRec', pred_time,...
              'delayTrans',pred_time,'MethodAssign',MethodAssign);
end
STE_shuffle = mean(tmp);

% Entropy of 'Target' signal to use as a normalization factor
H = quickEnt(symbRaster(1,:),'MethodAssign',{1,1,'Nat'});

%Define NSTE
NSTE = (STE - STE_shuffle)/H;
if NSTE < 0; NSTE = 0; end
end
