modpath = [cd '\LowerBody\MDHfiles\'];
c3dpath = [cd '\LowerBody\C3Dfiles\'];

standingdataname = [c3dpath 'Standing.c3d'];
templateFilename = [modpath 'Standing.mdh'];

classifierFilename = [datapath 'Classifier.mat'];

converttometers = 1;
reprocess = 1;
%Nearest neighbors to use as features (1 means only the marker of
%interest's position)
NumberOfNearestNeighbors = 1;
NumModels = 20; %Number of new models to generate

MaxNumRandomTargets = 3; %Adds randomness to number of tracking markers on each segment
jointnoise = 0.025; %Amount of noise (st. dev) for joint locations (meters)
 %Add some extra range to where the tracking markers can be as a ratio of
 %their actualy range.
extrarange = 0.1;
%Min distance a tracking marker can be from joint in vertical dimension
mindistfromjoint = 0.03;

%When training each classifier, enforce equality in number of positive
%instances of that class to negative instances.  The data set generally has
%way more negative than positive instances and seems to reduce performance
%of classifiers during cross-validation.
eliminateSampleBias = 1;

%Optimize SVM parameters using cross-validation error
SVMoptimization = 0;
%For each segment & joint label, this number chooses the number of random
%initial guesses to optimize over.  The best result is selected.
NumberOfInitialGuessesForSVMOptimization = 20;

if ~exist(classifierFilename,'file') || reprocess == 1;
    %% Get Standing data
    standingDataStruc = loadPhasespaceRecord(standingdataname, ...
        'c3dUnitsAreInMeters', 0,'forceReprocess',1);
    NumOfMarks = length(standingDataStruc.markerNames);
    fnames = fieldnames(standingDataStruc);
    for i = 1:length(fnames)
        if ~any(strcmp(fnames{i},{'markerNames' 'frameRate' 'time'}))
            d.(fnames{i}) = nanmean(standingDataStruc.(fnames{i}),1);
        end
    end
    
    x = double(cell2mat(struct2cell(d)));
    
    %% Get Labels
    [labels] = GetJointAndSegmentLabels(templateFilename);
    labs = cell2mat(struct2cell(labels));
    labelnames = fieldnames(labels);
    datanames = fieldnames(d);
    
    y = zeros(NumOfMarks,1);
    for i = 1:NumOfMarks
        y(i,1) = labs(strcmp(datanames{i},labelnames)==1);
    end
    
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
    models = x;
    xtrain = GetModelFeatures(x,NumberOfNearestNeighbors);
    ytrain = y;
    modeldex = [1;length(xtrain)];
    
    count = 0;
    for i = 1:NumModels-1
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
    
    %     rng(1e7); p = randperm(size(xtrain,1)); xtrain = xtrain(p,:);
    %     ytrain = ytrain(p);
    
    %% Training
    
    % Dimensionality reduction
    %     D = 12; [residuals,reconstructed] = pcares(xtrain,D); xtrain =
    %     reconstructed(:,1:D);
    
    %Train a unique classifier for each class
    Clabels = {'RightFoot' 'RightShank' 'RightThigh' 'LeftFoot' 'LeftShank' 'LeftThigh'...
        'Pelvis' 'RightAnkle' 'RightKnee' 'LeftAnkle' 'LeftKnee'};
    
    for i = 1:length(unique(ytrain))
        y = ytrain;
        x = xtrain;
        
        if eliminateSampleBias
            %Balance the number of labels
            NotTheClass = find(y~=i);
            dd = NotTheClass(1:end-length(find(y==i)));
            RemoveDex{i} = NotTheClass(dd);
            remlogical = ones(size(ytrain));
            remlogical(RemoveDex{i}) = 0;
            KeepDex{i} = find(remlogical==1);
            y(RemoveDex{i}) = [];
            x(RemoveDex{i},:) = [];
        end
        y(y~=i) = -1;
        y(y==i) = 1;
        yt{i} = y;
        xt{i} = x;
        
        fprintf(1,'Training %s classifier...\n',Clabels{i});
        Model = fitcsvm(xt{i},yt{i},'KernelFunction','rbf','Standardize',true);
        if SVMoptimization
            fprintf(1,'Cross validating to get Optimal Parameters...\n');
            
            %Set partition to be same each time
            c = cvpartition(length(yt{i}),'KFold',10);
            
            %Objective is to minimize cross-validated loss by tuning the
            %SVM paramters KernelFunction & Box Constraint
            minfn = @(z)kfoldLoss(fitcsvm(xt{i},yt{i},'CVPartition',c,...
                'KernelFunction','rbf','BoxConstraint',exp(z(2)),...
                'KernelScale',exp(z(1))));
            opts = optimset('TolX',5e-4,'TolFun',5e-4,'Display','final');
            
            
            %Perform 20 different optimizations using random IC and chosse
            %the best
            m = NumberOfInitialGuessesForSVMOptimization;
            fval = zeros(m,1);
            z = zeros(m,2);
            for j = 1:m;
                [searchmin, fval(j)] = patternsearch(minfn,randn(2,1),[],[],[],[],[],[],opts);
                z(j,:) = exp(searchmin);
            end
            z = z(fval == min(fval),:);
            
            %Final Model
            Model = fitcsvm(xt{i},yt{i},'KernelFunction','rbf',...
                'KernelScale',z(1),'BoxConstraint',z(2));
        end
        %     Classifiers.(Clabels{i}) = Model;
        fprintf(1,'Fitting posterior probabilities...\n');
        Classifiers.(Clabels{i}) = fitSVMPosterior(Model);
        
    end
    
    %% Testing
    labels = [];
    scores = [];
    for i = 1:length(unique(ytrain))
%             [L{i},S{i}] = predict(Classifiers.(Clabels{i}),xt{i});
        [labels(:,i),scores(:,:,i)] = predict(Classifiers.(Clabels{i}),xtrain);
    end
    
    
%     if ~eliminateSampleBias
%         scores = [];
%         labels = [];
%         for i = 1:length(L)
%             scores = cat(3,scores,S{i});
%             labels = cat(3,labels,L{i});
%         end
%         
%     else
%         for i = 1:length(unique(ytrain))
%             count = 0;
%             for j = 1:length(ytrain)
%                 if  ~any(j == RemoveDex{i})
%                     count = count+1;
%                     scores(j,:,i) = S{i}(count,:);
%                     labels(j,i) = L{i}(count);
%                 else
%                     scores(j,:,i) = [0 0];
%                     labels(j,i) = NaN;
%                 end
%             end
%         end
%     end
    
    posscores = scores;
    posscores(:,1,:) = [];
    posscores = squeeze(posscores);
    [Confidence,Y] = max(posscores,[],2);
    
    %% Save
    save(classifierFilename)
    
else
    load(classifierFilename);
end

wrongdex = find(Y~=ytrain);
fprintf(1,'Error Rate Identifying Markers = %.03g %%\n',sum(Y~=ytrain)/length(ytrain)*100);

if ~isempty(wrongdex)
    i = 1;
    
    figure
    h = find(modeldex == wrongdex(i));
    if isempty(h)
        ge = modeldex>wrongdex(i);
        h = find(ge==0,1,'last');
    end
    
    [~,coldex] = ind2sub(size(modeldex),h);
    x = models(modeldex(1,coldex):modeldex(2,coldex),:);
    y = ytrain(modeldex(1,coldex):modeldex(2,coldex));
    PlotMarkers(x,y);
    t = text(models(wrongdex(i),1),models(wrongdex(i),2),models(wrongdex(i),3),...
        sprintf('   Mislabeled %s\n   as %s',Clabels{ytrain(wrongdex(i))},Clabels{Y(wrongdex(i))}));
    set(t,'Color',[255 0 255]/255);
    
    plot3(models(wrongdex(i),1),models(wrongdex(i),2),models(wrongdex(i),3),...
        '.','MarkerSize',40,'Color',[255 0 255]/255)
end

