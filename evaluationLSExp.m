% Evaluation of OIS Forward Rate Curves
reloadData = false;
if ~exist('forwardCurvesLSExp','var') || reloadData
    load forwardCurveLSExp.mat
end

% settings
dt = 30/size(forwardCurvesLSExp,2);
lastPC = 4; %must be <= 9
% from this point and forward, the curve contains data for all maturities
% up to 30 years.

% The code block in "if reloadData" takes a long time to run. the outputs 
% have been saved in EVD...mat
if ~exist('totEigVal','var') && ~reloadData
    load EVD_LSExp.mat
end
if reloadData
    % contains all the forward rate curves
    forwardDiff = forwardCurvesLSExp(2:end,:) - forwardCurvesLSExp(1:end-1,:);
    forwardCov = cov(forwardDiff);
    [forwardEigVecLSExp, forwardEigValLSExp] = eig(forwardCov);
    forwardEigVecLSExp = fliplr(forwardEigVecLSExp);
    forwardEigValLSExp = rot90(forwardEigValLSExp,2);
    totEigValLSExp = sum(sum(forwardEigValLSExp));    
    loadingsPCLSExp = forwardEigVecLSExp(:,1:lastPC);
    weightsPCLSExp = diag(forwardEigValLSExp(1:lastPC,1:lastPC))/totEigValLSExp;
    totPCWeightLSExp = sum(weightsPCLSExp);
    save('EVD_LSExp.mat','loadingsPCLSExp','weightsPCLSExp','totPCWeightLSExp','totEigValLSExp');
end

plotPC = 4; % must be <= lastPC
legendPC = repmat(char(0),plotPC,3);
for i = 1:plotPC
    %tmpLeg = strcat('PC ',num2str(i),'(',num2str(round(weightsPC(i),4)*100),'%)');
    legendPC(i,:) = strcat('PC ',num2str(i));
end

figure(1);
hold on;
plot(dt:dt:size(loadingsPCLSExp,1)*dt,loadingsPCLSExp(:,1:plotPC));
legend(legendPC);