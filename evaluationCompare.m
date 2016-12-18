% Comparison between LSExp curves and LSExpSpline curves
reloadData = true;
compareForward = false;
comparePC = true;
invPC = 1; % when running with comparePC=true, if the factor loadings are 
           % inverted, just add or remove a minus sign from the 1 here.
dataLength = 30; % longest maturity (20 or 30)
var = 2; % set of spline curves

LSExpData = ['LSExp' num2str(dataLength) '.mat'];
splinesData = ['Splines' num2str(dataLength) '_' num2str(var) '.mat'];

forwardLSExpData = ['forwardCurve' LSExpData];
forwardSplinesData = ['forwardCurve' splinesData];
EVDLSExpData = ['EVD_' LSExpData];
EVDSplinesData = ['EVD_' splinesData];
if  reloadData
    load(forwardLSExpData);
    load(forwardSplinesData);
    load(EVDLSExpData);
    load(EVDSplinesData);
end

dtLSExp = dataLength/size(forwardCurvesLSExp,2);
dtSplines = dataLength/size(forwardCurvesSplines,2);

if compareForward
    figure(1);
    title('Forward Curves');
    for i = 1:size(forwardCurvesLSExp,1)
        plot(dtLSExp:dtLSExp:dataLength,forwardCurvesLSExp(i,:),...
            dtSplines:dtSplines:dataLength,forwardCurvesSplines(i,:));
        legend('LSExp','SplineApprox');
        pause(0.1);
    end
end

if comparePC
    figure(2);
    hold on;
    title('Loading Comparison: Shift');
    plot(dtLSExp:dtLSExp:dataLength,loadingsPCLSExp(:,1));
    plot(dtSplines:dtSplines:dataLength,invPC*loadingsPCSplines(:,1));
    msg1 = strcat('LSExp (',num2str(100*weightsPCLSExp(1)),'%)');
    msg2 = strcat('SplineApprox (',num2str(100*weightsPCSplines(1)),'%)');
    legend(msg1,msg2);
    
    figure(3);
    hold on;
    title('Loading Comparison: Twist');
    plot(dtLSExp:dtLSExp:dataLength,loadingsPCLSExp(:,2));
    plot(dtSplines:dtSplines:dataLength,invPC*loadingsPCSplines(:,2));
    msg1 = strcat('LSExp (',num2str(100*weightsPCLSExp(2)),'%)');
    msg2 = strcat('SplineApprox (',num2str(100*weightsPCSplines(2)),'%)');
    legend(msg1,msg2);
    
    figure(4);
    hold on;
    title('Loading Comparison: Butterfly');
    plot(dtLSExp:dtLSExp:dataLength,loadingsPCLSExp(:,3));
    plot(dtSplines:dtSplines:dataLength,invPC*loadingsPCSplines(:,3));
    msg1 = strcat('LSExp (',num2str(100*weightsPCLSExp(3)),'%)');
    msg2 = strcat('SplineApprox (',num2str(100*weightsPCSplines(3)),'%)');
    legend(msg1,msg2);
end
