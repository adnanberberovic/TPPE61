% Evaluation of OIS Forward Rate Curves
reloadData = false;
if ~exist('forwardMatrix','var') || reloadData
    load BlomvallLSexp1E0SEKTest11_cutoffEnd.mat
end

% settings
dt = 1/365;
lastPC = 4; %must be <= 9

allMaturities = 1;%2162+19; 
% from this point and forward, the curve contains data for all maturities
% up to 30 years.

% The code block in "if reloadData" takes a long time to run. the outputs 
% have been saved in IRS_Evaluation.mat
if ~exist('totEigVal','var')
    load OISf_cutoffEnd.mat
end
if reloadData
    % contains all the forward rate curves
    forwardMatrix20 = forwardMatrix(allMaturities:end,:);
    forwardDiff = forwardMatrix20(2:end,:) - forwardMatrix20(1:end-1,:);
    forwardCov = cov(forwardDiff);
    [forwardEigVec, forwardEigVal] = eig(forwardCov);
    forwardEigVec = fliplr(forwardEigVec);
    forwardEigVal = rot90(forwardEigVal,2);
    totEigVal = sum(sum(forwardEigVal));
    
end

loadingsPC = forwardEigVec(:,1:lastPC);
weightsPC = diag(forwardEigVal(1:lastPC,1:lastPC))/totEigVal;
totPCWeight = sum(weightsPC);

legendPC = repmat(char(0),lastPC,3);
for i = 1:lastPC
    %tmpLeg = strcat('PC ',num2str(i),'(',num2str(round(weightsPC(i),4)*100),'%)');
    legendPC(i,:) = strcat('PC ',num2str(i));
end

figure(1);
hold on;
plot(dt:dt:size(loadingsPC,1)*dt,loadingsPC);
legend(legendPC);