


dataAveragingWindow= 10;
smoothedLoad = movmean(load, 100);
Rsquared= [];
    for n = 50:10:1000
        fitEstimateWindow = n;
        base = smoothedLoad(1:fitEstimateWindow);
        [lineEstimate,S] = polyfit( extension(1:fitEstimateWindow), base,1);
        plotLineEstimate = lineEstimate(1)*extension+lineEstimate(2);
       %if lineEstimate(1)>0
        Rsquared(n) = 1 - (S.normr/norm(base - mean(base)))^2;
       %end
    end
    [bestRsquared, bestRsquaredIndex] = max(Rsquared);
    fitEstimateWindow = bestRsquaredIndex;
        base = smoothedLoad(1:fitEstimateWindow);
        [lineEstimate,S] = polyfit( extension(1:fitEstimateWindow), base,1);
        plotLineEstimate = lineEstimate(1)*extension+lineEstimate(2);
        initialLineFitRsquared = 1 - (S.normr/norm(base - mean(base)))^2
    for j=floor(dataAveragingWindow/2):(numel(extension)-dataAveragingWindow)
        testWindow = smoothedLoad(j-floor(dataAveragingWindow/2)+1:j+floor(dataAveragingWindow/2));
        lineWindow = plotLineEstimate(j-floor(dataAveragingWindow/2)+1:j+floor(dataAveragingWindow/2));
        %if ttest2(base, testWindow)
        if ttest2(lineWindow, testWindow)
        %if ttest(testWindow, plotLineEstimate(j))
        else
            divergencePoint = j;
        end
    end
    
    
    
    
    plot(extension, smoothedLoad);
    hold on;
    plot(extension,plotLineEstimate);
    scatter(extension(divergencePoint),plotLineEstimate(divergencePoint));
    hold off;
    
    pause(1);
    

