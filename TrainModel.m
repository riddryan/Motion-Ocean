N = 20; %Number of new models to generate from each seed model

%% Get Data
modpath = [cd '\LowerBody\MDHfiles\'];
c3dpath = [cd '\LowerBody\C3Dfiles\'];

mnames = {'Standing' 'Amy_normal' 'XiaoYu' ... 
           'Justin' 'Dana' 'Andrew' 'Alex'};
c3dnames = {'Standing' 'Stance' 'test1_combined' ...
             'test1_combined_3' 'test1_combined_4' 'test1_combined_5' ...
             'test2_combined'};
% mnames = {'Amy_normal'};
% c3dnames = {'Stance'};

modelnames = strcat(modpath,mnames,'.mdh');
c3dnames = strcat(c3dpath,c3dnames,'.c3d');


[xtrain,ytrain,models,modeldex] = GenerateTrainingData(c3dnames,modelnames,10);

%% Train the Classifier

C = TrainClassifier(xtrain,ytrain);

[Y,Confidence] = TestClassifier(C,xtrain,ytrain);

save([cd '\LowerBody\C_7seeds_10Noise_Opt.mat']);

