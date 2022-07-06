%% symbolicdata2states - converts data to states using ranking symbols
% Converts raw data to states using rankings of data values, similar to the
% symbolic entropy techniques applied in Staniek and Lehnertz, PRL 100,
% 2008, 158101. Note, wordstates will reduce the total number of time bins.
% Trial based data will be converted across trials at each individual time.
% Single trial data will be converted through time.
%
% Syntax: [StatesRaster] = symbolicdata2states(DataRaster, SymbInfo, varargin)
%
% Input:
%   DataRaster (cell array or double array) - trial data. If a double array
%       is used, it should be number of variables by number of time bins by
%       number of trials. If a cell array is used, it should have only one 
%       dimension and each element should be a double array with the 
%       dimensions listed above. Each element of the cell array is referred
%       to as a 'data category'.
%   SymbInfo (integer array) - sets the number of single time bin states to
%     combine into ranking sets. The array should have dimensions number of
%     data categories by 2. Each row controls the state conversion for a 
%     data category. The first column is the data category and the second
%     category is the number of contiguous single time bin states to
%     combine. If DataRaster contains only one data category, the first
%     column of SymbInfo can be removed. Note, all variables in a given
%     data category will be converted. Also, note that data categories not
%     listed in SymbInfo will not undergo state conversion.
%
% Variable Inputs:
%   (..., timeboundaries) - records the time boundaries for the time bins.
%     This must be a cell array or double array (matching DataRaster).
%     Each element (in the case of a cell array) should have dimension
%     number of time bins (matching the data category from DataRaster) by
%     2. 
%
%   (..., tau) - time delay variable (l in original paper; Staniek and 
%     Lehnertz, 2008). Number of samples between datapoints used for 
%     symbolization. Making this larger increases the time-scale of 
%     information transfer. Default is 1 (i.e. all neighboring points used). 
%
% Outputs:
%   StatesRaster (cell array or double array) - trial state data. If a
%     double array is used, it should be number of variables by number of
%     time bins by number of trials. Each element should be an integer
%     state number (state number = 1, 2, 3, ...). If a cell array is used,
%     it should have only one dimension and each element should be a double
%     array with the dimensions listed above. Each element of the cell
%     array is referred to as a 'data category'. StatesRaster will have
%     fewer time bins in the data categories converted via the symbolic
%     ranking proceedure.
%   timeboundaries (cell array or double array) - new time boundaries
%     following state conversion. If timeboundaries was not input, this
%     variable will be NaN. It will have the same structure as the input
%     timeboundaries, but it will have fewer bins for the data categories
%     that underwent state conversion.
%   
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See data2states for other stating options.
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% May 2016; Last revision: 13-May-2016


function [StatesRaster] = symbolicdata2states(DataRaster, SymbInfo, varargin)
%% Parse command line for parameters
TBOpt = false;
timeboundaries = NaN;
tau = 1;

iVarArg = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'timeboundaries',  timeboundaries = varargin{iVarArg+1}; iVarArg = iVarArg + 1; TBOpt = true;
        case 'tau',             tau = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(QUICKTE) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Perform Initial Operations

% Error check inputs and possible reorganize to ease later processing
if iscell(DataRaster)
    Back2ArrayFlag = false;
    if length(DataRaster) ~= 1
        if size(SymbInfo,2) ~= 2
            error('SymbInfo is not the appropriate size. Probably the data category is missing.')
        end
    else
        if size(SymbInfo,2) == 1
            SymbInfo = [ones(size(SymbInfo)),SymbInfo];
        elseif size(SymbInfo,2) ~= 2
            error('SymbInfo is not the appropriate size.')
        end
    end
    if TBOpt
        if length(timeboundaries) ~= length(DataRaster)
            error('timeboundaries and DataRasters do not have the same number of data categories.')
        else
            for iDC = 1:length(DataRaster)
                if size(timeboundaries{iDC},1) ~= size(DataRaster{iDC},2)
                    error(['timeboundaries and DataRaster are inconsistent in data category ',num2str(iDC)])
                end
            end
        end
    end
else
    Back2ArrayFlag = true;
    DataRaster = {DataRaster};
    if size(SymbInfo,2) == 1
        SymbInfo = [ones(size(SymbInfo)),SymbInfo];
    end
    if TBOpt
        if size(timeboundaries,1) ~= size(DataRaster{1},2)
            error('timeboundaries and DataRaster are inconsistent')
        end
    end
end

% Make sure we didn't get duplicate instructions
if length(unique(SymbInfo(:,1))) ~= length(SymbInfo(:,1))
    error('Duplicate data category instructions in SymbInfo.')
end



%% Perform the State Conversion

% Copy everything over
StatesRaster = DataRaster;

for iConv = 1:size(SymbInfo,1)
    
    if size(DataRaster{SymbInfo(iConv,1)},3) > 1 % Trial Based Data        
        % Combine states for unique rankings across trials
        nT = size(DataRaster{SymbInfo(iConv,1)},2);
        Temp1 = NaN(size(DataRaster{SymbInfo(iConv,1)}) - [0,SymbInfo(iConv,2) - 1,0]);
        for iVar = 1:size(DataRaster{SymbInfo(iConv,1)},1)
            for iT = 1:(nT - SymbInfo(iConv,2) + 1)
                Temp2 = squeeze(DataRaster{SymbInfo(iConv,1)}(iVar,iT:(iT + SymbInfo(iConv,2) - 1),:));
                [Y,Temp2] = sort(Temp2,1,'ascend');
                [B,I,Temp1(iVar,iT,:)] = unique(Temp2','rows');
            end
        end
        StatesRaster{SymbInfo(iConv,1)} = Temp1;
        
    else % Single Trial Data      
        % Combine states for unique rankings across time
        nT = size(DataRaster{SymbInfo(iConv,1)},2);
        Temp1 = NaN(size(DataRaster{SymbInfo(iConv,1)},1), ceil((nT - (SymbInfo(iConv,2)*tau)+1)/tau));
        for iVar = 1:size(DataRaster{SymbInfo(iConv,1)},1)
            Temp2 = NaN([size(Temp1,2),SymbInfo(iConv,2)]);
            for iT = 1:SymbInfo(iConv,2)
                Temp2(:,iT) = DataRaster{SymbInfo(iConv,1)}(iVar,iT*tau:tau:(nT - (SymbInfo(iConv,2)*tau)+iT*tau));
            end
            [~,ii]=sort(Temp2,2,'Ascend');
            [~,r]=sort(ii,2);
            [B,I,Temp1(iVar,:)] = unique(r,'rows');
        end
        StatesRaster{SymbInfo(iConv,1)} = Temp1;
        
    end
    
    if TBOpt % Correct the time boundaries
        
        timeboundaries{SymbInfo(iConv,1)}(1:(nT - SymbInfo(iConv,2) + 1),2) = timeboundaries{SymbInfo(iConv,1)}(SymbInfo(iConv,2):nT,2);
        timeboundaries{SymbInfo(iConv,1)}((nT - SymbInfo(iConv,2) + 2):nT,:) = [];
        
    end
    
end

% Convert StatesRaster back to an array, if necessary
if Back2ArrayFlag
    StatesRaster = StatesRaster{1};
end

