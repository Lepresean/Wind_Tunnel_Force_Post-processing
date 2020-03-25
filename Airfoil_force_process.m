%Post-processing for Wind Tunnel Airfoil Tests
%Sean Patrick Devey
%3_10_20
%% Description: top level prgm for postprocessing of data
%Purpose: reads in raw tab delimited data files generated by "LoadCell DAQ wForceCalc V4 SPD" LabVIEW prgm, converts to corrected aerodyanmic coefficients
%Requires Matlab functions: data_import.m, nondim_mom.m, nondim_force.m, wallCorr.m, coordTransform.m, Load2LD.m, horoBuoy.m, Mini45_Fx_Fy_error.m, nonDim.m, removeZeros.m
%files must not have .xls file extensions. Tab delimited .dat works well.
%reference data files are good to have as well, should be in .xlsx with columns: aoa, cx
%Units are assumed to be SI
clear; clc;
%% Initialization
chord = 8*2.54/100; %airfoil chord in meters
span = 11.8125*2.54/100; %airfoil span in meters new one is at 11+7/8, old at 11+13/16
num_pts = 10;

%here is where user imputs folder/list of files to post-process
current_folder = pwd;


% Microflaps runs --------------------------
% data_folder = strcat(current_folder,'\microflaps');
% output_folder = 'C:\Users\seanp\Documents\Box Sync\Everything\Documents\Research\Fall 2019 - Spring 2020\Data\post-processing\microflaps - no trip';
% theta = 151.6; %angle (deg) CCW offset mini45 x to airfoil

% Control runs -----------------------------
% data_folder = strcat(current_folder,'\control_mini45\11_18_19'); %control no trip
% data_folder = strcat(current_folder,'\control_mini45\1_21_20'); %control no trip
% output_folder = 'C:\Users\seanp\Documents\Box Sync\Everything\Documents\Research\Fall 2019 - Spring 2020\Data\post-processing\control - no trip';
% theta = 81.5;

data_folder = strcat(current_folder,'\control_mini45\1_25_20'); %control with trip
output_folder = 'C:\Users\seanp\Documents\Box Sync\Everything\Documents\Research\Fall 2019 - Spring 2020\Data\post-processing\control - trip';
theta = 81.5; %angle (deg) CCW offset mini45 x to airfoil
%% loop to do everything with each run
data_files = dir(fullfile(data_folder, '*.dat'));   %get list of each file in folder
% data_raw(:,:,1) = readmatrix(fullfile(data_folder,data_files(1).name),'OutputType','double');
for k = 1:length(data_files)    %read each file one at a time
    %% import data
    filename = fullfile(data_folder,data_files(k).name);
    data_raw = readmatrix(filename,'OutputType','double');
    %At this point, columns are assumed to be in this order:
    %time, speed (m/s), rho (kg/m3), AOA (deg), Fx (N)...Fz (N),	Tx (N-mm)...Tz(N-mm),
    %Load Cell 1 Mean Volt (V)...6 Mean Volt (V), Load Cell 1 Volt Std Dev (V)...6 Std Dev(V)
    
    data_raw(:,1,:) = [];   %trim off 1st column (timestamps)
    data_raw(:,10:21,:)=[]; %remove voltage info because I'm not going to use it anyway
    [row, col]=size(data_raw);
    
    %correct flow speed measurements at zero speed to zeros
    for i=1:row
        if data_raw(i,1)<2    %if flow speed is less than 2m/s
            data_raw(i,1)=0; %set flow speed to zero
        elseif isnan(data_raw(i,1)) %if flow speed is NaN (sometimes happens when tunnel is off)
            data_raw(i,1)=0; %set flow speed to zero
        end
    end
    
    %reduce data into mean and std by the groups at each data location
    row_new = row/num_pts;
    data = zeros(row_new,col);
    for j=1:col
        for n=1:row_new
            data(n,j) = mean(data_raw(num_pts*n-num_pts+1:n*num_pts,j));
            %std_data(n,j) = std(data(num_pts*n-num_pts+1:n*num_pts,j);
        end
    end
    %% loop for including bias correction tunnel corrections
    count=0;
    [~,name,~] = fileparts(filename);
    filename_dest = strcat(fullfile(output_folder,name),'.xlsx');
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    for biasCorr = [1,0]
        for tunnelCorr = [1,0]
            %% Convert to L,D
            data_processed = Load2LD(data,theta,chord,span,biasCorr,tunnelCorr);
            [row,col]=size(data_processed);
            %% save converted data to .dat files within subfolder
            if biasCorr == 0
                if tunnelCorr == 0
                    sheet_label = "no corrections";
                else
                    sheet_label = "tunnelCorr only";
                end
            else
                if tunnelCorr == 0
                    sheet_label = "biasCorr only";
                else
                    sheet_label = "biasCorr & tunnelCorr";
                end
            end
            writematrix(data_processed,filename_dest,'Sheet',sheet_label)
            %making sure to clear out any old data outside of the correct range
            writematrix(strings(25),filename_dest,'Sheet',sheet_label,'Range',strcat('A',num2str(row+1),':',alphabet(col),num2str(row+20)));
            writematrix(strings(25),filename_dest,'Sheet',sheet_label,'Range',strcat(alphabet(col+1),'1:',alphabet(col+10),num2str(row+20)));
        end
    end
end %the end of the loop for each single run read in from the specified folder
