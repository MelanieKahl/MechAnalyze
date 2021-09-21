%
% This program will determine the compressive modulus of a hydrogel from
% displacement - force data that are stored in an Excel file
% Peter Levett, 12 April, 2012

%% **************************************************************************
%%
clear all
tic

% Load the Excel file containing the displacement - force data. 
% The first sheet should contain the following:
%   Column A:   A list of gel numbers 
%   Column B:   The corresponding wet weight of each gel
%   Column C:   Gel heights, if known. For any gels for which height is not
%               known, this column should contain zero. 
%   Columns D onwards can contain notes to further identify the gels
%
% The following sheets should contain displacement - force data for all the
% gels listed in column A of the first sheet. 
%
% Stress - strain plots will be plotted for all gels, with a maxium of 4
% plots per figure. 

%%

disp('Welcome to Modfinder')


% Set the time and date for the filename
time = clock;
date = date;

hour = time(4);
minute = time(5);
if hour > 11
    append = 'pm';
else 
    append = 'am';
end
if minute <10
    minute = ['0' num2str(minute)];
end
hour = num2str(hour);
minute = num2str(minute);


file1 = ['Results', ' ' , date, ' ', hour,',', minute, append];
file2 = ['Results adjusted', ' ' , date, ' ', hour,',', minute, append];

% Ask the user to select an Excel file containing the data
[baseFileName, folder] = uigetfile('*.xl*', 'Please select your Excel file');
fullFileName = fullfile(folder, baseFileName);

% Load the information from the first sheet of the Excel file. 
gels = xlsread(fullFileName);
gel_numbers = gels(:, 1)';
gels_total = size(gel_numbers);
gels_total = gels_total(2);
Num_gels = gels_total;
Weights = gels(:, 2)';
Heights_input = gels(:, 3)';

% Create vectors of zeros to fill with heights, moduli and Rsq values. 
% A vector 'noise_end_vec' is also created, which will contain the indices
% of the final data points used to determine noise. This vector is not
% reported, but can be useful to know the regions where noise is
% calculated.
Modulus = zeros(1, Num_gels);
Rsq = zeros(1, Num_gels);
noise_end_vec = zeros(1, Num_gels);

full_figures = floor(Num_gels/4);
remainder = Num_gels - full_figures*4;

for i = 1:Num_gels
    load_gel_data = gel_numbers(i);
    num = num2str(load_gel_data);
    disp(['Loading data for gel number ' num])
    A{i} = xlsread(fullFileName, i + 1);
    i = i + 1;
end   

%% This section finds initial estimates for the heights and moduli. 

for i = 1:Num_gels
    i
% Load the weight and data for gel 'i'
    weight = Weights(i);
    data = A{i};
    legend_vals(i) = gel_numbers(i);
    %legend_vals_str = num2str(legend_vals);
    %i
    %x=input('stop');

% Separate the data into displacement and force vectors. Calculate the
% 'noise' from early data points, and 'zero' the force. Adjust the number
% of data points used to calculate noise based on where the gel starts.    
    disp = data(:,1);
    force = -data(:,2);
 
% Initially use data points 10 - 300 to calculate noise. This will be
% reduced as necessary. 
    noise_begin = 10;
    noise_end = 300;
    noise_test = 1;

    while noise_test < 2 
        noise = mean(force(noise_begin:noise_end)) ;
        %x=input('stop');
        force_tared = force-noise;

% Find the indices of the first 25 points that the force exceeds 0.0007 N.
% The height of the gel is defined as the first point that the tared force
% exceeds 0.0007 N four times consecutively. 
        start = find(force_tared > 0.0007, 200, 'first');
        diff_start = diff(start);
        counter = find(diff_start == 1, 1);
        counter1 = 1;

        while counter1 < 200
            if diff_start(counter + 1) == 1 && diff_start(counter + 2) == 1 && diff_start(counter + 3) ==1 && diff_start(counter + 4) ==1
                initial1 = start(counter + 1);
                height = disp(initial1);
                counter1 = 200;
            else
                counter = counter + 1;
            end
            counter1 = counter1 + 1;
        end
        if start(counter + 1) < (noise_end + 50)
            noise_end = noise_end - 10;
        else
            noise_test = 2;
        end
    end

% Over-ride the height found by the program with any height specified by
% the user in the Excel file.
    if Heights_input(i) ~= 0
        height = Heights_input(i);
        initial1 = find(disp > height, 1, 'last');
        height = disp (initial1);
    end

% Store the found/entered height into the vector 'Heights' for reporting.
    Heights(i) = height;

% Store the index of the entry used as the final point to calculate noise.
    noise_end_vec(i)= noise_end;

% Display the following error message if the program is using data points
% from 'within' the gel to calculate noise.
    if start(counter + 1) < noise_end
        x = input('Error: Please reduce the number of points for noise');
    end
%--------------------------------------------------------------------------

% Calculate the strain and stress vectors. Determine the modulus as the 
% slope of the curve beteween 10 and 15% strain.

    disp_new = disp(initial1:end);
    force_new = force_tared(initial1:end);

    strain = (height - disp_new)/height;
    strain_store{i} = strain;
    
    area = weight/height;
    stress = force_new/area;
    stress_store{i} = stress;
    initial_strain_index = find(strain < 0.1, 1, 'last');
    final_strain_index = find(strain > 0.15, 1, 'first');
    strain_region = strain(initial_strain_index:final_strain_index);
    stress_region = stress(initial_strain_index:final_strain_index);


% Use polyfit to fit a linear model to the data of stress and strain 
% between 10 - 15% strain. Determine the modulus as the slope of the 
% line of best fit.
    p = polyfit(strain_region, stress_region, 1);
    Modulus_kPa = p(1)*1000;
    Modulus(i) = Modulus_kPa;
    
% Use polyval to determine the output of the model that has been fitted,
% and corrcoef to determine the correlation between the data and the model.
    output = polyval(p, strain_region);
    Correlation = corrcoef(stress_region, output);
    Rsq(i) = Correlation(1, 2)*Correlation(1,2);
end
%%

%% Produce stress-strain plots for each gel
% Plot data for a maxiumum of four graphs on each set of axes. 
c = 1;
f = 1;
while f <= full_figures
    figure
    plot(strain_store{c},stress_store{c},'b', strain_store{c+1},stress_store{c+1},'k',strain_store{c+2},stress_store{c+2},'r',strain_store{c+3},stress_store{c+3},'g')
    legend(['Gel ' num2str(legend_vals(c))] ,[ 'Gel ' num2str(legend_vals(c+1))], ['Gel ' num2str(legend_vals(c+2))],[ 'Gel ' num2str(legend_vals(c+3))] )
    xlabel('Strain')
    ylabel('Stress (kPa)')
    c = c+4;
    f = f+1;
end
if remainder == 3
    figure
    plot(strain_store{i-2},stress_store{i-2},'b',strain_store{i-1},stress_store{i-1},'k',strain_store{i},stress_store{i},'r')
    legend(['Gel ' num2str(legend_vals(i-2))], ['Gel ' num2str(legend_vals(i-1))],[ 'Gel ' num2str(legend_vals(i))] )
    xlabel('Strain')
    ylabel('Stress (kPa)')
elseif remainder == 2
    figure
    plot(strain_store{i-1},stress_store{i-1},'b',strain_store{i},stress_store{i},'k')
    legend(['Gel ' num2str(legend_vals(i-1))],[ 'Gel ' num2str(legend_vals(i))] )
    xlabel('Strain')
    ylabel('Stress (kPa)')
elseif remainder == 1
    figure
    plot(strain_store{i},stress_store{i},'b')
    legend([ 'Gel ' num2str(legend_vals(i))] )
    xlabel('Strain')
    ylabel('Stress (kPa)')
end   
%--------------------------------------------------------------------------

%% Produce a table containing the estimated data
figure
% Convert the row vectors to column vectors for displaying in a table
gel_numbers_col = gel_numbers';
Modulus_col = Modulus';
Weights_col = Weights';
Heights_col = Heights';
Rsq_col = Rsq';
m = [ gel_numbers_col Weights_col Heights_col Modulus_col Rsq_col ];

cnames = {'Gel Number','Weight (mg)', 'Height (mm)', 'Modulus (kPa)', 'R sq'};
m_lab = {cnames; m};
t = uitable('Data', m, 'columnName', cnames, 'columnWidth', {'auto'}, 'Position', [100 20 410 300]);

%%
% Check if the user is happy with the outputs, or if any gels should be
% measured again. 
%
%%
overall_check = 1;
overall_check2 = 1;
while overall_check < 2
    check = input('Would you like to test any gels again? 1 = yes, 2 = no :') ;
    while overall_check2 < 2
        if check == 1 
            num_check = input('Please enter the gel number of the gel you would like to check: ');
            num_check_index = find(num_check == gel_numbers);
            
            weight_loop = Weights(num_check_index);
            data_loop = A{num_check_index};
            disp_loop = data_loop(:,1);
            force_loop = -data_loop(:,2);
            
            
            % Initially use data points 10 - 300 to calculate noise. This will be
            % reduced as necessary. 
            noise_begin_loop = 10;
            noise_end_loop = 300;
            noise_test_loop = 1;

            while noise_test_loop < 2 
                noise_loop = mean(force_loop(noise_begin_loop:noise_end_loop)) ;
                force_tared_loop = force_loop-noise_loop;

            % Find the indices of the first 25 points that the force exceeds 0.0007 N.
            % The height of the gel is defined as the first point that the tared force
            % exceeds 0.0007 N four times consecutively. 
                start_loop = find(force_tared_loop > 0.0007, 200, 'first');
                diff_start_loop = diff(start_loop);
                counter = find(diff_start_loop == 1, 1);
                counter1 = 1;

                while counter1 < 200
                    if diff_start_loop(counter + 1) == 1 && diff_start_loop(counter + 2) == 1 && diff_start_loop(counter + 3) ==1 && diff_start_loop(counter + 4) ==1
                        initial1_loop = start_loop(counter + 1);
                        height_loop = disp_loop(initial1_loop);
                        counter1 = 200;
                    else
                        counter = counter + 1;
                    end
                    counter1 = counter1 + 1;
                end
                if start_loop(counter + 1) < (noise_end_loop + 50)
                    noise_end_loop = noise_end_loop - 10;
                else
                    noise_test_loop = 2;
                end
            end
             
% Produce a plot of displacement vs force to allow the user to judge if the
% value for height that has been determined is reasonable. Plot the graph
% between 'begin_disp' and 'end_disp'. 
                if initial1_loop > 200
                    begin_disp_loop = initial1_loop - 200;
                elseif initial1_loop > 150
                    begin_disp_loop = initial1_loop - 150;
                elseif initial1_loop > 100
                    begin_disp_loop = initial1_loop - 100  ;
                else 
                    begin_disp_loop = initial1_loop - 20;
                end
                end_disp_loop = initial1_loop + 250;
                strain_plot_loop = disp_loop(begin_disp_loop : end_disp_loop);
                stress_plot_loop = force_tared_loop(begin_disp_loop : end_disp_loop);
                stress_min_loop = min(stress_plot_loop);
                stress_max_loop = max(stress_plot_loop);
                figure
                plot(strain_plot_loop, stress_plot_loop, [height_loop height_loop], [stress_min_loop-.0015 stress_max_loop+.001])
                title('Use this graph to determine the height of the gel')
                xlabel('Displacement')
                ylabel('Force')
                height_check = height_loop

                d = input('Does the height in the graph look ok? 1 = yes, 2 = no: ');
                    if d == 2
                        height_loop = input('Please enter the gel height (mm): ');
                        initial1_loop = find(disp_loop > height_loop, 1, 'last');
                        height_loop = disp_loop(initial1_loop);
                    end

                disp_new_loop = disp_loop(initial1_loop:end);
                force_new_loop = force_tared_loop(initial1_loop:end);

                strain_loop = (height_loop - disp_new_loop)/height_loop;
                area_loop = weight_loop/height_loop;
                stress_loop = force_new_loop/area_loop;
                initial_strain_index_loop = find(strain_loop < 0.1, 1, 'last');
                final_strain_index_loop = find(strain_loop > 0.15, 1, 'first');
                strain_region_loop = strain_loop(initial_strain_index_loop:final_strain_index_loop);
                stress_region_loop = stress_loop(initial_strain_index_loop:final_strain_index_loop);

                Heights(num_check_index) = height_loop;

                p_loop = polyfit(strain_region_loop, stress_region_loop, 1);

                Modulus_kPa = p_loop(1)*1000
                Modulus(num_check_index) = Modulus_kPa;

                output_loop = polyval(p, strain_region_loop);
                Correlation_loop = corrcoef(stress_region_loop, output_loop);
                Rsq(num_check_index) = Correlation_loop(1, 2)*Correlation_loop(1,2);
                

            final_check = input('Would you like to test another gel? 1 = yes, 2 = no: ') ;            
            if final_check == 2
                overall_check2 = 2;                
            end
        else
            
            overall_check2 = 2;
            overall_check = 2;
        end
    overall_check = 2;
    end
end

% Determine if any values have changed, and if so, produce an updated
% table. Export initial estimates and final results to Excel. 
m2 = [ gel_numbers' Weights' Heights' Modulus' Rsq' ];
minimum_val = min(m./m2);
minimum_val = min(minimum_val);
maximum_val = max(m./m2);
maximum_val = max(maximum_val);
if maximum_val > 1.001
    minimum_val = 0.5;
elseif maximum_val < 0.999
    minimum_val = 0.5;
end
if minimum_val < 0.999
    minimum_val = 0.5;
end

if minimum_val < 0.999  
    col_width = size(cnames);
    col_width = col_width(2);
    m3 = zeros(1, col_width);
    m4 = [m; m3; m2];
       
    figure
    t = uitable('Data', m4, 'columnName', cnames, 'columnWidth', {'auto'}, 'Position', [100 20 410 300]);
    xlswrite(file2, m4)
else
    xlswrite(file1,m);
end
toc
%% The End