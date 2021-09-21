%getCompressionModulus v3.3.7 (differentiation mode)
%extracting compressive modulus and sample height from piston compression
%v3.3 adding more flexible raw data selection method
%v3.3.1 added more flexible method to specify crossectional area
%v3.3.2 general GUI cleanup
%v3.3.3 config file access debugging
%v3.3.4 added adaptive raw data detection and adaptive slope thresholding
%and added a check if the provided force and position vectors have the
%correct direction with an option to invert them
%v3.3.5 reverted to non-adaptive detection method
%v3.3.6 added manual specification of sample height
%v3.3.7 added exponential fit output 

%Dominik Schneidereit
%Institute of Medical Biotechnology
%13.07.2020

disp('MechAnalyze v3.3.7');
disp('Dominik Schneidereit');
disp('Institute of Medical Biotechnology');

%start of code
clear;
%get home folder
pathEnvironment = mfilename('fullpath');

%initialize config file with default values if it does not exist
%otherwhise load the config data
global filterString;
global lastLoadPath;
global lastSavePath;
configFileName = [char(java.lang.System.getProperty('user.home')) '\MechAnalyzeConfig.csv'];
if not(isfile(configFileName))
    configData = table();
    configData(1,1) = {'C:\'};
    configData.Properties.VariableNames(1) = {'lastSelectedPath'};
    configData(1,2) = {''};
    configData.Properties.VariableNames(2) = {'lastFilterString'};
    configData(1,3) = {'C:\'};
    configData.Properties.VariableNames(3) = {'lastSavePath'};
    configData(1,4) = {'C:\'};
    configData.Properties.VariableNames(4) = {'lastLoadPath'};
else
    configData = readtable(configFileName,'Delimiter',',');
end
%due to bugs i decided to remove the function that remembers the filter
%string for now
%filterString = configData.lastFilterString{:};
filterString = '';
lastLoadPath = configData.lastLoadPath{:};
lastSavePath = configData.lastSavePath{:};
%file handling
basePath = uigetdir(configData.lastSelectedPath{:}, 'Select search-basepath');
if basePath == 0
    return
end
global areaCrossecValues;
global specHeightValues;
global fileList;
%fileList = getfn(basePath,'Specimen_RawData');
fileList = getfn(basePath,'.csv|.xls|.xlsx');
global filteredFileList;
filteredFileList = fileList(1);
global bottomOffset;
bottomOffset = 50;
global bOK
bOK = 0;

%make the user filter for and select the data sets for evaluation
if not(isempty(filterString))
    initiallist = fileList(find(contains(fileList, filterString)));
    filteredFileList = initiallist;
else
    initiallist = fileList;
end
f = figure;
f.NumberTitle = 'off';
f.Name = 'Select data set(s) to evaluate';
f.ToolBar = 'none';
f.MenuBar = 'none';
oldPos = f.Position;
f.Position = [oldPos(1),oldPos(2),700,420];
f.SizeChangedFcn = @adaptTableSize_callback;
% maybe better use listdlg?
uil = uicontrol(f, 'Style', 'listbox');
uil.String = initiallist;
uil.Max = length(fileList);
uil.Position = [0 bottomOffset f.Position(3) f.Position(4)-bottomOffset];
uil.Callback = @listSelection_callback;

filtField = uicontrol(f,'Style','edit','Callback', @filtButton_callback, 'String',filterString);
filtField.Position = [filtField.Position(1)+60 filtField.Position(2) filtField.Position(3)+100 filtField.Position(4)];
filtTag = uicontrol(f,'Style','text', 'Position', [filtField.Position(1)-80 filtField.Position(2)-3 80 filtField.Position(4)], 'String', 'Filter String:');
bfilt = uicontrol(f,'Style','pushbutton','Callback',@filtButton_callback);
bfilt.String = 'Apply';
bfilt.Position = [filtField.Position(1)+160 filtField.Position(2) 60 filtField.Position(4)];
bcfilt = uicontrol(f,'Style','pushbutton','Callback',@clearFiltButton_callback);
bcfilt.String = 'Clear';
bcfilt.Position = [filtField.Position(1)+220 filtField.Position(2) 60 filtField.Position(4)];
bc = uicontrol(f,'Style','pushbutton','Callback',@closebutton_callback);
bc.Position =  [filtField.Position(1)+520 filtField.Position(2) 60 filtField.Position(4)];
bc.String = 'OK';

drawnow;
waitfor(f);

if bOK == 0
    disp('Process aborted');
    return
end

fileList = filteredFileList;

global bButtonPressed;
bButtonPressed = 0;
%global selected_cells;
%selected_cells = [];
global figureHandle
figureHandle = [];
global extensionInitPos;
global loadInitPos;
extensionInitPos = [];
loadInitPos = [];

%get data environment
parentFolderNames = cell(numel(fileList),1);
fileNames = cell(numel(fileList),1);
filePathList = cell(numel(fileList),1);
fileIndexList = cell(numel(fileList),1);
dataSetNames = cell(numel(fileList),1);

%get file environment
disp([num2str(numel(fileList)) ' Files selected for processing']);
for i=1:numel(fileList)
    [folderPath, fileName, fileExtension] = fileparts(fileList{i});
    fileNames(i,1) = {fileName};
    splitFolderPath = split(folderPath,'\');
    parentFolderNames(i,1) = splitFolderPath(end);
    filePathList(i,1) = {folderPath};
end

%specify a UID for all data
uniqueId = join([parentFolderNames fileNames], ';');

%check for potential ambiguities in the file structure
if not(numel(unique(uniqueId))==numel(uniqueId))
    questdlg('Some datasets seem to have the same combination of file and folder name, please rename them or analyze them seperately','Error','OK','OK');
    return
end

%populate results array
resultsArray = [parentFolderNames fileNames];

%define crossection mode
crossecMode = questdlg('How would you like to specify the sample crossection area?','Select crossection mode','Constant for all data sets','Specify individually for each data set','Cancel','Constant for all data sets');

%obtain crossection data, using either one user specified value or reading
%values from a user specified excel sheet
switch crossecMode
    case 'Cancel'
        return
    case 'Constant for all data sets'
        rawAreaCrossec = inputdlg('Please specify sample crossectional area (in mm²)','Specify Area',1);    
        areaCrossec = str2double(rawAreaCrossec);
        if (isempty(areaCrossec))
            disp('You did not enter a number (e.g. 8.314)');
            disp('Ending program run..');
            return
        elseif (isnan(areaCrossec))
            disp('You did not enter a number (e.g. 8.314)');
            disp('Ending program run..');
            return
        else
            resultsArray(:,3) = {areaCrossec};
        end    
    case 'Specify individually for each data set'
        
        bOK = 0;
        f = figure;
        f.NumberTitle = 'off';
        f.Name = 'Enter the crossection data or load/save from file';
        f.ToolBar = 'none';
        f.MenuBar = 'none';
        f.SizeChangedFcn = @adaptTableSize_callback;
        uit = uitable(f);
        uit.Position = [0 bottomOffset f.Position(3) f.Position(4)-bottomOffset];
        areaCrossecValues = num2cell(zeros(size(fileNames,1),size(fileNames,2)));
        uit.Data = [parentFolderNames fileNames areaCrossecValues];
        uit.ColumnEditable = [false false true];
        uit.ColumnName = {'Folder Name'; 'File Name'; 'Crossection area (mm²)'};
        uit.ColumnWidth = [{(length(parentFolderNames{1})+3)*6} {(length(fileNames{1})+3)*6}];
        bLoad = uicontrol(f,'Style','pushbutton','Callback',@loadButton_callback);
        bLoad.String = 'Load from File';
        bLoad.Position = [bLoad.Position(1)+0 bLoad.Position(2) 80 bLoad.Position(4)];
        bSave =uicontrol(f,'Style','pushbutton','Callback',@saveButton_callback);
        bSave.String = 'Save to File';
        bSave.Position = [bSave.Position(1)+100 bSave.Position(2) 80 bSave.Position(4)];
        bOKbutton = uicontrol(f,'Style','pushbutton','Callback',@okCrossec_callback);
        bOKbutton.String = 'OK';
        bOKbutton.Position = [bOKbutton.Position(1)+360 bOKbutton.Position(2) 60 bOKbutton.Position(4)];
        drawnow;
        waitfor(f);
        
        if bOK == 0
            disp('Process aborted');
            return
        end
        
        resultsArray(:,3) = areaCrossecValues;                   
end

%query height mode
heightMode = questdlg('Would you like to autodetect the height of all samples or specify the heigh of some samples manually?','Select height mode','Autodetect all','Specify heights manually','Cancel','Autodetect all');

%populate height value array if desired
autoHeight = 0;
switch heightMode
    case 'Cancel'
        return
    case 'Autodetect all'
        autoHeight = 1;
    case 'Specify heights manually'
        autoHeight = 0;
        bOK = 0;
        f = figure;
        f.NumberTitle = 'off';
        f.Name = 'Enter the crossection data or load/save from file';
        f.ToolBar = 'none';
        f.MenuBar = 'none';
        f.SizeChangedFcn = @adaptTableSize_callback;
        uit = uitable(f);
        uit.Position = [0 bottomOffset f.Position(3) f.Position(4)-bottomOffset];
        specHeightValues = num2cell(zeros(size(fileNames,1),size(fileNames,2)));
        uit.Data = [parentFolderNames fileNames specHeightValues];
        uit.ColumnEditable = [false false true];
        uit.ColumnName = {'Folder Name'; 'File Name'; 'Sample height (mm)'};
        uit.ColumnWidth = [{(length(parentFolderNames{1})+3)*6} {(length(fileNames{1})+3)*6}];
        bLoad = uicontrol(f,'Style','pushbutton','Callback',@loadButton_callback);
        bLoad.String = 'Load from File';
        bLoad.Position = [bLoad.Position(1)+0 bLoad.Position(2) 80 bLoad.Position(4)];
        bSave =uicontrol(f,'Style','pushbutton','Callback',@saveButton_callback);
        bSave.String = 'Save to File';
        bSave.Position = [bSave.Position(1)+100 bSave.Position(2) 80 bSave.Position(4)];
        bOKbutton = uicontrol(f,'Style','pushbutton','Callback',@okSpecHeight_callback);
        bOKbutton.String = 'OK';
        bOKbutton.Position = [bOKbutton.Position(1)+360 bOKbutton.Position(2) 60 bOKbutton.Position(4)];
        drawnow;
        waitfor(f);
        
        if bOK == 0
            disp('Process aborted');
            return
        end
end

%display first data file and make user select correct columns
f = figure;
f.NumberTitle = 'off';
f.Name = 'Specify the start of load/extension data';
f.ToolBar = 'none';
f.MenuBar = 'none';
f.SizeChangedFcn = @adaptTableSize_callback;
uit = uitable(f);
uit.Position = [0 bottomOffset f.Position(3) f.Position(4)-bottomOffset];
origDisplayData = table2cell(readtable(fileList{1}));
uit.Data = origDisplayData;
global iselect;
iselect = 0;
uit.CellSelectionCallback = @cellselection_callback;
bc = uicontrol(f,'Style','pushbutton','Callback',@closebutton_callback);
bc.Position = [100 20 60 20];
bc.String = 'OK';
bOK = 0;
bc.Enable = 'off';
br = uicontrol(f,'Style','pushbutton','Callback',@redobutton_callback);
br.String = 'Redo';
infoField = uicontrol(f, 'Style', 'pushbutton', 'Enable', 'inactive', 'String', '<html><center> Select the first <b><font color="blue">extension</font></b> value, then select the first <b><font color="green">load</font></b> value in the data set and confirm with <b>OK</b></center></html>');
infoField.Position = [250 10 300 35];

drawnow;

waitfor(f);
if bOK == 0
    disp('Process aborted');
    return
end

%update config file
configData.lastLoadPath = {lastLoadPath};
configData.lastSavePath = {lastSavePath};
configData.lastSelectedPath = {basePath};
configData.lastFilterString = {filterString};
writetable(configData,configFileName,'Delimiter',',');

bFlipExtension = false;
bflipLoad = false;

%iterate over each dataset and perform evaluation routines
f = figure;
for i= 1:numel(fileList)
    %read data from files
    dataTable = readtable(fileList{i});
    %if data is not numeric, try to convert it, otherwise use the numeric
    %values that are not flagged as NaN
    if not(isnumeric(dataTable{extensionInitPos(1),extensionInitPos(2)}))
        extensionNaN = str2double(dataTable{extensionInitPos(1):end,extensionInitPos(2)});
        extension = extensionNaN(not(isnan(extensionNaN)));
        loadNaN = str2double(dataTable{loadInitPos(1):end,loadInitPos(2)});
        load = loadNaN(not(isnan(loadNaN)));
    else
        extension = dataTable{extensionInitPos(1):end,extensionInitPos(2)}(not(isnan(dataTable{extensionInitPos(1):end,extensionInitPos(2)})));
        load = dataTable{loadInitPos(1):end,loadInitPos(2)}(not(isnan(dataTable{loadInitPos(1):end,loadInitPos(2)})));
    end
    
    %flip data if requested
    if i == 1
        if mean(extension)>0
            decisionToFlip = questdlg(['The detected extension values seem to be positive.' newline 'However, negative extension values are expected with coordinate 0' newline 'being the collision position of piston and surface.' newline 'Do you want to flip the coordinate system? (multiply values with -1)'],'Warning!','Yes','No','Cancel','Yes');
            switch decisionToFlip
                case 'Yes'
                    bFlipExtension = true;
                case 'No'
                    bFlipExtension = false;
                case 'Cancel'
                    return;
            end
        end
        
        if mean(load)<0
            decisionToFlip = questdlg(['The detected load values seem to be mainly negative.' newline 'However, compressive force load is expected to be positive.' newline 'Do you want to flip the coordinate system? (multiply values with -1)'],'Warning!','Yes','No','Cancel','Yes');
            switch decisionToFlip
                case 'Yes'
                    bflipLoad = true;
                case 'No'
                    bflipLoad = false;
                case 'Cancel'
                    return
            end
        end
    end
    
    if (bFlipExtension)
        extension = -extension;
    end
    
    if (bflipLoad)
        load = -load;
    end
    
    %filter low level data noise
    smoothLoad = movmean(load(1:end-10),20);
    
    %get relevant points from dataset
    [maxLoad, iMaxLoad] = max(load(1:end-10));
    [minLoad, iMinLoad] = min(smoothLoad(1:iMaxLoad));
    
    %calculate load threshold
    deltaLoad = maxLoad-minLoad;
    loadThreshold = 0.0005*deltaLoad;
    
    %estimate divergence point for sample height determination
    lookupWindow = 20;   
    slopes = diff(smoothLoad)./diff(extension(1:end-10));
    smoothSlopes = movmean(slopes,lookupWindow);
    %to allow broader experimental variability, the slope threshold is
    %scaled to the data set
    %divergencePoint = find(smoothSlopes(lookupWindow:end)>(0.02), 1)+lookupWindow;
    %scaledSlopethreshold = 0.02 * (max(smoothSlopes)-smoothSlopes(2));
    scaledSlopethreshold = 0.02;
    
    %if height is specified by the user, use that data to determine
    %divergence point
    if isequal(specHeightValues(i),{0})||isnan(specHeightValues{i})
        divergencePoint = find(smoothSlopes(lookupWindow:end)>(scaledSlopethreshold), 1)+lookupWindow;
    else
        divergencePoint = find(extension>-specHeightValues{i},1);
    end
    %display plot to illustrate the divergence point
    f.NumberTitle = 'off';
    f.Name = [parentFolderNames{i} '\' fileNames{i}];
    %f.ToolBar = 'none';
    %f.MenuBar = 'none';
    sp1 = subplot(2,2,1);
    plot(extension, load);
    hold on;
    scatter(extension(divergencePoint),load(divergencePoint));
    hold off;
    axis([(min(extension)-0.5) (max(extension)+0.5) (minLoad-loadThreshold) (((maxLoad-minLoad-loadThreshold)/4)+(minLoad-loadThreshold))])
    %axis([(min(extension)-0.5) (max(extension)+0.5) (minLoad-loadThreshold) maxLoad])
    xlabel('Extension (mm)')
    ylabel('Load (N)')

    %calculate sample height from divergence point
    L0 = extension(divergencePoint);
    sampleHeight = -L0;
    F0 = load(divergencePoint);
    
    %calculate strain and stress
    strain = (extension-L0)./sampleHeight;
    
    %check for a crossection area entry for this sample, if none is
    %present, use the average of all specified values and flag a warning
    if not(isempty(resultsArray{i,3}))
        crossecArea=resultsArray{i,3};
        bCrossecIsSpecified = 1;
    else
        crossecArea=nanmean(cell2mat(resultsArray(:,3)));
        disp('Warning! Missing Area crossection information, using average of specified values..');
        bCrossecIsSpecified = 0;
    end   
    stress = (load-F0)*1000/crossecArea;
    strainStress = [strain stress];
    
    %check for failure 
    failureCase = 0;
    %slopes = diff(smoothLoad)./diff(extension);
    smoothSlopes = movmean(slopes(1:end-2),lookupWindow);
    firstSlopeDip = find(smoothSlopes<-(0.05*(max(smoothSlopes)-smoothSlopes(2))), 1);
    if not(isempty(firstSlopeDip))
        failureCase = 1;
        [failureLoad,indexShift] = max(load((firstSlopeDip-lookupWindow):(firstSlopeDip+lookupWindow)));
        iFailureLoad = firstSlopeDip-lookupWindow+indexShift;
        failureStress = stress(iFailureLoad);
        failureStrain = strain(iFailureLoad);
    end
        
    %determine compressive modulus by courve fitting to strain stress
    lowerLimit = 0.1;
    upperLimit = 0.15;
    relevantRegion = strainStress((lowerLimit<strainStress(:,1)&strainStress(:,1)<upperLimit), :);
    compModulusFit = polyfit( relevantRegion(:,1), relevantRegion(:,2),1);
    compModFitPlot = compModulusFit(1)*strain+compModulusFit(2);
    compressiveModulus = compModulusFit(1); %kPa
    
    %determine exponential compressive modulus by courve fitting to strain stress
    
    if failureCase==1
        extendedRelevantRegion = strainStress((0<strainStress(:,1)&strainStress(:,1)<failureStrain), :);
    else
        extendedRelevantRegion = strainStress((0<strainStress(:,1)&strainStress(:,1)<strain(end)), :);
    end
    [expCompModFitCourve, expCompModFitGoF] = fit(extendedRelevantRegion(:,1), extendedRelevantRegion(:,2),'exp1');
    expCompModConfInt = confint(expCompModFitCourve);
    expCompModStdDev = (expCompModConfInt(2,:)-expCompModConfInt(1,:))./4;
    expCompMod = [expCompModFitCourve.a expCompModFitCourve.b];
    %log results
    resultsArray(i,4) = {compressiveModulus};
    resultsArray(i,5) = {sampleHeight};
    resultsArray(i,8) = {expCompModFitCourve.a};
    resultsArray(i,9) = {expCompModStdDev(1)};
    resultsArray(i,10) = {expCompModFitCourve.b};
    resultsArray(i,11) = {expCompModStdDev(2)};
    resultsArray(i,12) = {expCompModFitGoF.rsquare};
    
    %display stress strain courve with fit courve
    sp2 = subplot(2,2,2);
    plot(strain,stress,'Color','black')
    hold on;
    plot(strain, expCompModFitCourve.a*exp(strain.*expCompModFitCourve.b),'Color','blue');
    plot(strain, compModFitPlot,'Color','red')
    hold off;
    xlabel('Strain (-)')
    ylabel('Stress (kPa)')
    CurrentAxes.XLim(1) = 0;
    CurrentAxes.YLim = [min(stress) max(stress)];
    
    %display detected failure case
    sp3 = subplot(2,2,3);
    cla(sp3);
    axis on;
    switch failureCase
        case 1
            plot(extension,load)
            hold on;
            scatter(extension(iFailureLoad), failureLoad)
            hold off;
            failureText = 'Failure stress: ';
            failureText1 = 'Failure strain: ';
            failureDisplayStress = failureStress;
            failureDisplayStrain = failureStrain;       
        case 0
            axis off;
            displayText1 = {'No mechanical failure'};
            if exist('textItem1')
                delete(textItem1);
            end
            textItem1 = text(0,0.5,displayText1);
            failureLoad = NaN;
            failureStress = NaN;
            failureStrain = NaN;
            failureText = 'No failure until stress: ';
            failureText1 = 'No failure until strain: ';
            failureDisplayStress = stress(end);
            failureDisplayStrain = strain(end); 
    end
    xl1 = xlabel('Extension (mm)');
    yl1 = ylabel('Load (N)');
    resultsArray(i,6) = {failureStress};
    resultsArray(i,7) = {failureStrain};
    
    %display results in the graph
    sp4 = subplot(2,2,4);
    cla(sp4);
    axis off;
    if bCrossecIsSpecified == 1
        displayText = {['CompressiveModulus: ' num2str(compressiveModulus,'% .2f') ' kPa'], ['Sample height: ' num2str(sampleHeight,'% .2f') ' mm'], [failureText num2str(failureDisplayStress,'% .2f') ' kPa'],[failureText1 num2str(failureDisplayStrain,'% .3f') ]};
    else
        displayText = {['CompressiveModulus: ' num2str(compressiveModulus,'% .2f') ' kPa'], ['Sample height: ' num2str(sampleHeight,'% .2f') ' mm'], [failureText num2str(failureDisplayStress,'% .2f') ' kPa'],[failureText1 num2str(failureDisplayStrain,'% .3f') ], 'Warning! No area value was specified for this sample!','Average of specified values is used!'};
    end
    if exist('textItem')
        delete(textItem);
    end
    textItem = text(0,0.5,displayText);
    disp(['File ' num2str(i) ' of ' num2str(numel(fileList)) ' || Sample height: ' num2str(sampleHeight,'% .2f') 'mm  Compressive modulus: ' num2str(compressiveModulus,'% .2f') 'kPa'])
    drawnow;
    % saving the figures as PNG and matlab FIG
    print([filePathList{i} '\' fileNames{i} 'Eval.png'],'-dpng','-r400');
    savefig([filePathList{i} '\' fileNames{i} 'Eval.fig']);
    % !! this is how long the resulting plots are displayed set this to 0
    % for fastest evaluation run and to disable updating plots, set to higher numers
    % (in seconds) get more time to look at the results
    % pause(0.0001);
    
end

%save results to file
resultsTable = cell2table(resultsArray);
resultsTable.Properties.VariableNames = [{'FolderName'} {'FileName'} {'CrossecArea_mmSquared'} {'compressiveModulus_kPa'} {'SampleHeight_mm'}, {'FailureStress_kPa'},{'FailureStrain'},{'ExpCompressiveModulus_Amplitude_kPa'},{'ExpCompressiveModulus_Amplitude_StdDev_kPa'},{'ExpCompressiveModulus_Exponent'},{'ExpCompressiveModulus_Exponent_StdDev'},{'ExpCompressiveModulus_FitRsquare'}];
writetable(resultsTable, [basePath '\EvalResults.csv'], 'Delimiter',',');

%message at end of run with data structure details
disp('Evaluation is done.');
disp([num2str(numel(fileList)) ' data files were processed.']);
disp('Find the results in "EvalResults.csv".');
disp(['in the folder ' basePath]);
disp('The figures are stored alongside the raw data files.');

questdlg(['Evaluation routine completed' newline 'Find the results in "' basePath '\EvalResults.csv".'],'Done','OK','OK');

close(f);


function okSpecHeight_callback(src,event)
    global bOK;
    global specHeightValues;
    bOK = 1;
    
    specHeightValues = src.Parent.Children(end).Data(:,3);
    close(src.Parent);
end


function okCrossec_callback(src,event)
    global bOK;
    global areaCrossecValues;
    bOK = 1;
    
    areaCrossecValues = src.Parent.Children(end).Data(:,3);
    close(src.Parent);
end

function loadButton_callback(src, event)
    global lastLoadPath;
    [fileName, filePath, idx] = uigetfile([lastLoadPath '*.xlsx']);
    if numel(fileName) == numel(filePath) == numel(idx) && idx == 0
        disp('Loading cancelled');
    else
        loadedData = table2cell(readtable([filePath fileName],'ReadVariableNames',false));
        loadedDataNoempty = [];
        
        %check if the loaded data has at least 3 columns        
        if size(loadedData,2)<3
            questdlg(['Incompatible data format.' newline 'Please make sure the specified spreadsheet contains data pertaining to files selected for evaluation and it is in the is in the expected format:' newline '(folder name|file name|crossection area)'],'Error','OK','OK');
            return
        end
        
        %remove incomplete rows from loaded data
        for i=1:size(loadedData,1)
            if not(isempty(loadedData{i,1}))&&not(isempty(loadedData{i,2}))&&not(isnan(loadedData{i,3}))
                loadedDataNoempty = cat(1,loadedDataNoempty, loadedData(i,:));
            end
        end
        loadedData = loadedDataNoempty;
        %bContainsRelevantData = sum(ismember(loadedData(:,1),src.Parent.Children(end).Data(:,1)),'all')
        %bContainsRelevantData2 = sum(ismember(loadedData(:,2),src.Parent.Children(end).Data(:,2)),'all')
        if or(sum(ismember(loadedData(:,1),src.Parent.Children(end).Data(:,1)),'all')==0, sum(ismember(loadedData(:,2),src.Parent.Children(end).Data(:,2)),'all')==0)
            questdlg(['No relevant data could be fount in the specified file.' newline 'Please make sure the specified spreadsheet contains data pertaining to files selected for evaluation and it is in the is in the expected format:' newline '(folder name|file name|crossection area)'],'Error','OK','OK');
        else
            loadedUID = join([loadedData(:,1) loadedData(:,2)],'');
            originalUID = join([src.Parent.Children(end).Data(:,1) src.Parent.Children(end).Data(:,2)],'');
            if numel(loadedUID)==numel(unique(loadedUID))
                for i=1:numel(originalUID)
                    for j=1:numel(loadedUID)
                        if isequal(originalUID(i),loadedUID(j))
                            src.Parent.Children(end).Data(i,3) = loadedData(j,3);
                        end

                    end
                end
                lastLoadPath = filePath;
            else
                %disp('There are duplicate entries in the specified file, please remove duplicates or select a different file. Loading was aborted');
                questdlg('There are duplicate entries in the specified file, please remove duplicates or select a different file. Loading was aborted','Error','OK','OK');
            end
        end
        
            
    end
end

function saveButton_callback(src, event)
    global lastSavePath;
    [fileName, filePath, idx] = uiputfile([lastSavePath '*.xlsx']);
    if numel(fileName) == numel(filePath) ==numel(idx) && idx==0
        disp('Saving cancelled');
    else
        writetable(cell2table(src.Parent.Children(end).Data),[filePath fileName], 'WriteVariableNames', false);
        lastSavePath = filePath;
    end
end



function listSelection_callback(src, event)
    global filteredFileList
    filteredFileList = src.String(src.Value);
end

function clearFiltButton_callback(src,event)
    global fileList; 
    src.Parent.Children(end).String = fileList;
    src.Parent.Children(end-1).String = '';
end

function filtButton_callback(src,event)
    global fileList;
    global filterString;
    filterString = src.Parent.Children(end-1).String;
    src.Parent.Children(end).String = fileList(find(contains(fileList, filterString)));
    src.Parent.Children(end).Value = 1;
end


function cellselection_callback(src,event)
    %global selected_cells
    global iselect
    global extensionInitPos
    global loadInitPos
    iselect = iselect+1;
    if not(isempty(event.Indices))
        selected_cells = event.Indices;
        if iselect == 1
            highlightcolor = '"Blue"';
            extensionInitPos = selected_cells;
        else
            highlightcolor = '"Green"';
            loadInitPos = selected_cells;
        end
        if not(isnumeric(src.Data{selected_cells(2),selected_cells(2)}))
            for i=selected_cells(1):size(src.Data,1)
                src.Data(i,selected_cells(2)) ={['<html><body font color=' highlightcolor '>' src.Data{i,selected_cells(2)} '</body></html>']};
            end
        else
            for i=selected_cells(1):size(src.Data,1)
                %newContent = ['<html><body font color=' highlightcolor '>' num2str(src.Data{i,selected_cells(2)}) '</body></html>']
                src.Data(i,selected_cells(2)) ={['<html><body font color=' highlightcolor '>' num2str(src.Data{i,selected_cells(2)}) '</body></html>']};
            end
        end
            
        
    end
    if iselect > 3
        src.Parent.Children(end-1).Enable = 'on';
        src.Enable = 'inactive';
    end
%     uitjf = findjobj_fast(src);
%     uitjfv = uitjf.getViewport.getView;
%     uitjfv.changeSelection(size(src.Data,1)-1,selected_cells(1)-1,false,false);
%     uitjfv.changeSelection(selected_cells(2)-1,selected_cells(1)-1,true,true);
end





function redobutton_callback(src,event)
    %global origDisplayData;    
    %src.Parent.Children(end).Data = origDisplayData;
    global fileList;
    src.Parent.Children(end).Data = table2cell(readtable(fileList{1}));
    src.Parent.Children(end).Enable = 'on';
    src.Parent.Children(end-1).Enable = 'off';
    global iselect;
    iselect = 0;
end

function closebutton_callback(src,event)
    close(src.Parent)
    global bOK
    bOK = 1;
end

function adaptTableSize_callback(src,event)
    global bottomOffset;
    src.Children(end).Position = [0 bottomOffset src.Position(3) src.Position(4)-bottomOffset];
end
