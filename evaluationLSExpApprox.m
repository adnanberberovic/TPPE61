% Evaluation of OIS Forward Rate Curves with LSExpSpline
reloadData = true;
dataLength = 30; % longest maturity in years
var = 1; % set of spline curves
pvar = 100; % penalty induced
dataFile = ['Splines' num2str(dataLength) '_' num2str(var) '_' num2str(pvar) '.mat'];

forwardData = ['forwardCurve' dataFile];
pcaData = ['EVD_' dataFile];
if ~exist('forwardCurvesSplines','var') || reloadData
    load(forwardData);
end

% settings
dt = dataLength/size(forwardCurvesSplines,2);
lastPC = 4; %must be <= 9
% from this point and forward, the curve contains data for all maturities
% up to 30 years.

% The code block in "if reloadData" takes a long time to run. the outputs 
% have been saved in EVD...mat
if ~exist('totEigVal','var') && ~reloadData
    load(pcaData);
end
if reloadData
    % contains all the forward rate curves
    forwardDiff = forwardCurvesSplines(2:end,:) - forwardCurvesSplines(1:end-1,:);
    forwardCov = cov(forwardDiff);
    [forwardEigVecSplines, forwardEigValSplines] = eig(forwardCov);
    forwardEigVecSplines = fliplr(forwardEigVecSplines);
    forwardEigValSplines = rot90(forwardEigValSplines,2);
    totEigValSplines = sum(sum(forwardEigValSplines));
    loadingsPCSplines = forwardEigVecSplines(:,1:lastPC);
    weightsPCSplines = diag(forwardEigValSplines(1:lastPC,1:lastPC))/totEigValSplines;
    totPCWeightSplines = sum(weightsPCSplines);  
    save(pcaData,'loadingsPCSplines','weightsPCSplines','totPCWeightSplines','totEigValSplines');
end

plotPC = 4; % must be <= lastPC
legendPC = repmat(char(0),plotPC,3);
for i = 1:plotPC
    legendPC(i,:) = strcat('PC ',num2str(i));
end

figure(1);
hold on;
plot(dt:dt:size(loadingsPCSplines,1)*dt,loadingsPCSplines(:,1:plotPC));
legend(legendPC);