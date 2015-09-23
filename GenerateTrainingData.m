function [xtrain,ytrain,models,modeldex] = GenerateTrainingData(c3dnames,modelnames,N,varargin)
% GENERATETRAININGDATA generates labeled features associated with human
% skeletal models.
% 
% [xtrain,ytrain] = GENERATETRAININGDATA(c3dnames,modelnames) returns
% features XTRAIN and labels for those features YTRAIN using Visual3D model
% files MODELNAMES built from corresponding .c3d motion files C3DNAMES
% 
% [xtrain,ytrain,models,modeldex] =
% GENERATETRAININGDATA(c3dnames,modelnames,N) returns MODELS which have the
% raw model marker data generated and MODELDEX which specifies which
% indices of models correspond to single markers.  N specifies the number
% of random models to randomly generate from each seed model in MODELNAMES.
%
%
if length(c3dnames)~=length(modelnames)
    error('Must have equal Number of .c3d & .mdh files')
end

if nargin == 2
    N = 20;
end

%Nearest neighbors to use as features (1 means only the marker of
%interest's position)
NumberOfNearestNeighbors = 1;
MaxNumRandomTargets = 3; %Adds randomness to number of tracking markers on each segment
jointnoise = 0.025; %Amount of noise (st. dev) for joint locations (meters)
 %Add some extra range to where the tracking markers can be as a ratio of
 %their actualy range.
extrarange = 0.1;
%Min distance a tracking marker can be from joint in vertical dimension
mindistfromjoint = 0.03;


opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'jointnoise'
        jointnoise = val;
      case 'mindistfromjoint'
          mindistfromjoint = val;
      case 'extrarange'
          extrarange = val;
      case 'MaxNumRandomTargets'
          MaxNumRandomTargets = val;
      case 'NumberOfNearestNeighbors'
          NumberOfNearestNeighbors = val;
      
      
    otherwise
      error('\nError: incorrect varargin\n')
  end
end


Mnum = length(c3dnames);

xtrain = [];
ytrain = [];
models = [];
modeldex = [];
for m = 1:Mnum
    %% Get Standing data
    S = loadPhasespaceRecord(c3dnames{m},'c3dUnitsAreInMeters',0,'forceReprocess',1);
    
    fnames = fieldnames(S);
    d = [];
    for i = 1:length(fnames)
        if ~any(strcmp(fnames{i},{'markerNames' 'frameRate' 'time' 'analogs'}))
            d.(fnames{i}) = nanmean(S.(fnames{i}),1);
        end
    end

    %% Get Labels
    [labels] = GetJointAndSegmentLabels(modelnames{m});
    labs = cell2mat(struct2cell(labels));
    
    y = labs;
    labelnames = fieldnames(labels);
    NumOfMarks = length(labelnames);
    
    datanames = fieldnames(d);
    for i = 1:length(datanames)
        if ~any(strcmp(datanames{i},labelnames))
            d = rmfield(d,datanames{i});
        end
    end
    
    dordered = [];
    for i = 1:length(labelnames)
       dordered.(labelnames{i}) = d.(labelnames{i}); 
    end
    d = dordered;
%     
    x = double(cell2mat(struct2cell(d)));
    
    rng(1000);
    P = randperm(length(y));
    y = y(P);
    x = x(P,:);
    
    xorig = x;
    yorig = y;
    
    %% Generate training data
    % For the seed model, generate new models by 1) adding a random number
    % of tracking markers for each segment at uniformally random locations
    % between each joint. 2) New joint markers that have normally
    % distributed noise from their original locations.
    %
    % The scale and rotation of the model is taken care of in the feature
    % extraction function GetModelFeatures.
    models = [models; x];
    xtrain = [xtrain; GetModelFeatures(x,NumberOfNearestNeighbors)];
    ytrain = [ytrain; y];
    modeldex = [modeldex [1;length(xtrain)]];
    
    for i = 1:N-1
        x = xorig;
        y = yorig;
        
        %% Noisy Joint Locations
        %         noise =
        %         bsxfun(@minus,bsxfun(@times,jointnoise,rand(2,3)),jointnoise/2);
        rng(i);
        noise = jointnoise*randn(2,3);
        %Right ankle
        x(y==8,:) = x(y==8,:) + noise;
        
        %Right Knee
        x(y==9,:) = x(y==9,:) + noise;
        
        %Left ankle
        x(y==10,:) = x(y==10,:) + noise;
        
        %Left knee
        x(y==11,:) = x(y==11,:) + noise;
        
        
        %% Random Tracking Targets
        rightankmax = max(x(y==8,:));
        rightkneemax = max(x(y==9,:));
        leftankmax = max(x(y==10,:));
        leftkneemax = max(x(y==11,:));
        rightankmin = min(x(y==8,:));
        rightkneemin = min(x(y==9,:));
        leftankmin = min(x(y==10,:));
        leftkneemin= min(x(y==11,:));
        
        %Right Foot
        rng(i+5000);
        T = ceil(rand(1,1)*MaxNumRandomTargets)+1;
        pmin = min(x(y==1 | y==8,:),[],1);
        pmax = max(x(y==1 | y==8,:),[],1);
        prange = range(x(y==1 | y==8,:),1);
        b = pmax + extrarange*prange/2;
        a = pmin - extrarange*prange/2;
        b(3) = rightankmin(3); %;Don't allow feet markers above ankle
        randmarkers = bsxfun(@plus,a,bsxfun(@times,b-a,rand(T,3)));
        x(y==1,:) = []; y(y==1) = [];
        x = [x; randmarkers];
        y = [y; 1*ones(T,1)];
        
        %Right Shank
        rng(i+5000*2);
        T = ceil(rand*MaxNumRandomTargets);
        pmin = min(x(y==2 | y==8 | y==9,:),[],1);
        pmax = max(x(y==2 | y==8 | y==9,:),[],1);
        prange = range(x(y==2 | y==8 | y==9,:),1);
        b = pmax + extrarange*prange/2;
        a = pmin - extrarange*prange/2;
        b(3) = rightkneemin(3) - mindistfromjoint; %Don't allow shank markers above knee
        a(3) = rightankmax(3) + mindistfromjoint; %Don't allow shank markers below ankle
        randmarkers = bsxfun(@plus,a,bsxfun(@times,b-a,rand(T,3)));
        x(y==2,:) = []; y(y==2) = [];
        x = [x; randmarkers];
        y = [y; 2*ones(T,1)];
        
        %Right Thigh
        rng(i+5000*3);
        T = ceil(rand*MaxNumRandomTargets)+1;
        pmin = min(x(y==3 | y==9,:),[],1);
        pmax = max(x(y==3 | y==9,:),[],1);
        prange = range(x(y==3 | y==9,:),1);
        b = pmax + extrarange*prange/2;
        a = pmin - extrarange*prange/2;
        a(3) = rightkneemax(3) + mindistfromjoint; %Don't allow thigh markers below knee
        randmarkers = bsxfun(@plus,a,bsxfun(@times,b-a,rand(T,3)));
        x(y==3,:) = []; y(y==3) = [];
        x = [x; randmarkers];
        y = [y; 3*ones(T,1)];
        
        %Left Foot
        rng(i+5000*4);
        T = ceil(rand*MaxNumRandomTargets)+1;
        pmin = min(x(y==4 | y == 10,:),[],1);
        pmax = max(x(y==4 | y == 10,:),[],1);
        prange = range(x(y==4 | y == 10,:),1);
        b = pmax + extrarange*prange/2;
        a = pmin - extrarange*prange/2;
        b(3) = rightankmin(3); %Don't allow feet markers above ankle
        randmarkers = bsxfun(@plus,a,bsxfun(@times,b-a,rand(T,3)));
        x(y==4,:) = []; y(y==4) = [];
        x = [x; randmarkers];
        y = [y; 4*ones(T,1)];
        
        %Left Shank
        rng(i+5000*5);
        T = ceil(rand*MaxNumRandomTargets);
        pmin = min(x(y==5 | y == 10 | y == 11,:),[],1);
        pmax = max(x(y==5 | y == 10 | y == 11,:),[],1);
        prange = range(x(y==5 | y == 10 | y == 11,:),1);
        b = pmax + extrarange*prange/2;
        a = pmin - extrarange*prange/2;
        b(3) = leftkneemin(3) - mindistfromjoint; %Don't allow shank markers above knee
        a(3) = leftankmax(3) + mindistfromjoint; %Don't allow shank markers below ankle
        randmarkers = bsxfun(@plus,a,bsxfun(@times,b-a,rand(T,3)));
        x(y==5,:) = []; y(y==5) = [];
        x = [x; randmarkers];
        y = [y; 5*ones(T,1)];
        
        %Left Thigh
        rng(i+5000*6);
        T = ceil(rand*MaxNumRandomTargets)+1;
        pmin = min(x(y==6 | y == 11,:),[],1);
        pmax = max(x(y==6 | y == 11,:),[],1);
        prange = range(x(y==6 | y == 11,:),1);
        b = pmax + extrarange*prange/2;
        a = pmin - extrarange*prange/2;
        a(3) = leftkneemax(3) + mindistfromjoint; %Don't allow thigh markers below knee
        randmarkers = bsxfun(@plus,a,bsxfun(@times,b-a,rand(T,3)));
        x(y==6,:) = []; y(y==6) = [];
        x = [x; randmarkers];
        y = [y; 6*ones(T,1)];
        
        %Right Pelvis
        rng(i+5000*7);
        T = ceil(rand*MaxNumRandomTargets) + 2;
        pmin = min(x(y==7,:),[],1);
        pmax = max(x(y==7,:),[],1);
        prange = range(x(y==7,:),1);
        b = pmax + extrarange*prange/2;
        a = pmin - extrarange*prange/2;
        randmarkers = bsxfun(@plus,a,bsxfun(@times,b-a,rand(T,3)));
        x(y==7,:) = []; y(y==7) = [];
        x = [x; randmarkers];
        y = [y; 7*ones(T,1)];
        
        [~,I] = sort(y);
        x = x(I,:);
        y = y(I,:);
        
        xtrain = [xtrain; GetModelFeatures(x,NumberOfNearestNeighbors)];
        ytrain = [ytrain; y];
        models = [models; x];
        modeldex = [modeldex [modeldex(end)+1 ;modeldex(end) + length(y)]];
    end
    
end


end