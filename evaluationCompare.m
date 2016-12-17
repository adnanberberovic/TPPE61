% Comparison between LSExp curves and Spline curves
reloadData = true;
compareForward = false;
comparePC = true;

if  reloadData
    load forwardCurveLSExp.mat
    load forwardCurveSplines.mat
    load EVD_LSExp.mat
    load EVD_Splines.mat
end

dtLSExp = 30/size(forwardCurvesLSExp,2);
dtSplines = 30/size(forwardCurvesSplines,2);

if compareForward
    figure(1);
    title('Forward Curves');
    for i = 1:size(forwardCurvesLSExp,1)
        plot(dtLSExp:dtLSExp:30,forwardCurvesLSExp(i,:),...
            dtSplines:dtSplines:30,forwardCurvesSplines(i,:));
        legend('LSExp','SplineApprox');
        pause(0.1);
    end
end

if comparePC
    figure(2);
    hold on;
    title('Loading Comparison: Shift');
    plot(dtLSExp:dtLSExp:30,loadingsPCLSExp(:,1));
    plot(dtSplines:dtSplines:30,-loadingsPCSplines(:,1));
    msg1 = strcat('LSExp (',num2str(100*weightsPCLSExp(1)),'%)');
    msg2 = strcat('SplineApprox (',num2str(100*weightsPCSplines(1)),'%)');
    legend(msg1,msg2);
    
    figure(3);
    hold on;
    title('Loading Comparison: Twist');
    plot(dtLSExp:dtLSExp:30,loadingsPCLSExp(:,2));
    plot(dtSplines:dtSplines:30,-loadingsPCSplines(:,2));
    msg1 = strcat('LSExp (',num2str(100*weightsPCLSExp(2)),'%)');
    msg2 = strcat('SplineApprox (',num2str(100*weightsPCSplines(2)),'%)');
    legend(msg1,msg2);
    
    figure(4);
    hold on;
    title('Loading Comparison: Butterfly');
    plot(dtLSExp:dtLSExp:30,loadingsPCLSExp(:,3));
    plot(dtSplines:dtSplines:30,-loadingsPCSplines(:,3));
    msg1 = strcat('LSExp (',num2str(100*weightsPCLSExp(3)),'%)');
    msg2 = strcat('SplineApprox (',num2str(100*weightsPCSplines(3)),'%)');
    legend(msg1,msg2);
end
