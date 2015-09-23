function [C] = TrainClassifier(xtrain,ytrain)
%When training each classifier, enforce equality in number of positive
%instances of that class to negative instances.  The data set generally has
%way more negative than positive instances and seems to reduce performance
%of classifiers during cross-validation.
eliminateSampleBias = 1;

%Optimize SVM parameters using cross-validation error
SVMoptimization = 1;
%For each segment & joint label, this number chooses the number of random
%initial guesses to optimize over.  The best result is selected.
NumberOfInitialGuessesForSVMOptimization = 5;

Clabels = {'RightFoot' 'RightShank' 'RightThigh' 'LeftFoot' 'LeftShank' 'LeftThigh'...
    'Pelvis' 'RightAnkle' 'RightKnee' 'LeftAnkle' 'LeftKnee'};

%% Train a binary SVM Classifier for each element of Clabels
for i = 1:length(unique(ytrain))
    y = ytrain;
    x = xtrain;
    
    %% Eliminate Sample Bias
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
    
    %% Fit Classifier
    fprintf(1,'Training %s classifier...\n',Clabels{i});
    Model = fitcsvm(x,y,'KernelFunction','rbf','Standardize',true);
    %% Optimize Parameters
    if SVMoptimization
        fprintf(1,'Cross validating to get Optimal Parameters...\n');
        
        %Set partition to be same each time
        c = cvpartition(length(y),'KFold',10);
        
        %Objective is to minimize cross-validated loss by tuning the SVM
        %paramters KernelFunction & Box Constraint
        minfn = @(z)kfoldLoss(fitcsvm(x,y,'CVPartition',c,...
            'KernelFunction','rbf','BoxConstraint',exp(z(2)),...
            'KernelScale',exp(z(1))));
        opts = optimset('TolX',5e-4,'TolFun',5e-4,'Display','final');
        
        
        %Perform different optimizations using random IC and choose the
        %best
        m = NumberOfInitialGuessesForSVMOptimization;
        fval = zeros(m,1);
        z = zeros(m,2);
        for j = 1:m;
            [searchmin, fval(j)] = patternsearch(minfn,randn(2,1),[],[],[],[],[],[],opts);
            z(j,:) = exp(searchmin);
        end
        z = z(fval == min(fval),:);
        
        %Final Model
        Model = fitcsvm(x,y,'KernelFunction','rbf',...
            'KernelScale',z(1),'BoxConstraint',z(2));
    end
    
    %% Get Posterior Probabilities
    fprintf(1,'Fitting posterior probabilities...\n');
    C.(Clabels{i}) = fitSVMPosterior(Model);
    
end

%% Get Test Error

end