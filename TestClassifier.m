function [Y,Confidence] = TestClassifier(C,x,y)
if nargin==2
    y = [];
end

labels=[]; scores = [];
names = fieldnames(C);
for i = 1:length(names)
    [labels(:,i),scores(:,:,i)] = predict(C.(names{i}),x);
end

posscores = scores;
posscores(:,1,:) = [];
posscores = squeeze(posscores);
[Confidence,Y] = max(posscores,[],2);

if ~isempty(y)
    fprintf(1,'Error Rate Identifying Markers = %.03g %%\n',sum(Y~=y)/length(y)*100);
end

end