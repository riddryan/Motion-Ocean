function features = GetModelFeatures(markers,K)
if nargin == 1
K = 5; %Number of nearest neighbors to use
end

% Get position of each marker relative to center of position
meanpos = mean(markers,1);
markers = bsxfun(@minus,markers,meanpos);

%Normalize scale of data
posnorm = zeros(size(markers,1),1);
for i = 1:size(markers,1)
    posnorm(i) = norm(markers(i,:));
end
scalefactor = 1/mean(posnorm);
markers = scalefactor * markers;

%Rotate the markers such that they have the least amount of covariance in
%the x & y directions. (Make the x & y directions independent)
PCAVectors = pca(markers(:,1:2));
Theta = -atan2(PCAVectors(2,1),PCAVectors(1,1));
Rotation = [cos(Theta) -sin(Theta); sin(Theta) cos(Theta)];
markers(:,1:2) = (Rotation*markers(:,1:2)')';

%Find nearest neighbors
IDX = knnsearch(markers,markers,'K',K);

% For each marker, calculate ratio of other markers that have greater
% values in the x, y, and z dimensions.
NumOfMarks = length(markers);
RelativeRatios = zeros(size(markers));
for i = 1:NumOfMarks
    greatherthandex = bsxfun(@ge,markers,markers(i,:));
    RelativeRatios(i,:) = (sum(greatherthandex,1) - 1)/(NumOfMarks - 1);
    % MAYBE ALSO GET MEAN POSITION OF ALL MARKERS ABOVE, RIGHT, AND FORWARD
    % OF CURRENT MARKERfd
end

%Prepare features for classification
features = [reshape(markers(IDX,:),NumOfMarks,K*3) RelativeRatios];


end