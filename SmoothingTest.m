

for i = 1:1:20
    
    smoothedLoad = movmean(logLoad, i);
    base = smoothedLoad(1:dataAveragingWindow);
    lineEstimate = polyfit( extension(1:dataAveragingWindow), base,1);
    plotLineEstimate = lineEstimate(1)*extension+lineEstimate(2);
    
    plot(extension, smoothedLoad);
    hold on;
    plot(extension,plotLineEstimate);
    hold off;
    pause(1);
    i
end
