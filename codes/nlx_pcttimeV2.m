clear all
% Paths and user defined inputs
addpath(genpath('/Users/aliizadi/documents/work/codes/theta analysis tools/new codes/'))
data_dir = '/Users/aliizadi/desktop/analysis2/';
analyze_times_directory =  '/Users/aliizadi/documents/work/codes/timesaves/';

% EEG Properties
plot_low_freq = 6;
plot_high_freq = 10;
freqs = 6:10;
wavelet_cycles = 6;
duration = 3;
srate = 1000;

Automatic_Start_Input = false;    %Allows the EEG start times to be entered automaticly (Line: )
%end user defined inputs

[files] = dir([data_dir,'*.ncs']);
file_total = length(files);

%% % ------------------------- Time Stamps ------------------------- % %%
%Obtain Time Frames
fprintf('\nObtaining Start and End Times\n');
End_Thresh = 10;                          %The amount of time allowed to truncate before asking for a different start and end time
ID_Column  = 'B';                         %Column letter with the subject IDs
Row_Limit = '200';                        %The number of rows it will check before it stops

if Automatic_Start_Input
    id_loc = [ID_Column,'1:',ID_Column,Row_Limit];
    [ndata,text,alldata] = xlsread(Start_Time_File,Start_Time_Sheet,id_loc);
    id_array = alldata;
    total_rows = length(id_array);
    
    start_loc = [Start_Time_Column,'1:',Start_Time_Column,Row_Limit];
    [ndata,text,alldata] = xlsread(Start_Time_File,Start_Time_Sheet,start_loc);
    start_array = alldata;
end

for I = 1:file_total
    %Load File
    file_to_analyze = files(I).name;
    short_fname = strtok(file_to_analyze,'.'); %set up timestamp fname
    %%%%%% end_temp (end time of EEG)
    
    %Get ID (finds the Rat's ID)

    
    %gets start and end time
    if Automatic_Start_Input
        %Check if ID is on the Excel sheet and get row
        row_loc = 0;
        full_pass = true;
        msgId = 'MATLAB:nonIntegerTruncatedInConversionToChar';
        warnStruct = warning('off',msgId);
        for J = 1:total_rows
            if length(char(id_array{J}))~= length(char(fID))
                continue
            elseif char(id_array{J}) == char(fID);
                row_loc = J;
                full_pass = false;
                break
            end
        end
        warning(warnStruct);
        if full_pass
            check=inputdlg({sprintf('There is no subject %s in the Excel file\nClick OK to skip or Cancel to exit.',fID)},'Subject Check',1,{'1'});
            if check{1} ~= 1
            end
            continue
        end      
        start_time = start_array{row_loc};
        end_time = (start_time+Test_Dur);
        %Automatic EEG start input 
        if end_time > end_temp
            if End_Thresh < (end_time-end_temp)
                while true
                    warning('The EEG Ended before the recorded end time at %fs\nManually put in Start and End times for %s',end_temp,file_to_analyze);
                    answer = inputdlg({sprintf('%s\nselect start time (s)',file_to_analyze);'select end time (s)'},file_to_analyze,1,{'1';'2'});
                    start_time = str2num(answer{1});
                    end_time =   str2num(answer{2});
                    if end_time <= end_temp
                        break
                    end
                    warning('Those times are still invalid.');
                end
            else
                end_time = end_temp;
                warning('The EEG Ended before the recorded end time truncating %fs',(endtime-end_temp));
            end
        end
    else
        %Manual EEG start input
        cd(analyze_times_directory)
        if exist([short_fname,'_timesave.mat']) ~=2
            answer = inputdlg({sprintf('%s\nselect start time (s)',file_to_analyze);...
                 'select end time (s)'},file_to_analyze,1,{'1';'2'});
            start_time = str2num(answer{1});
            end_time =   str2num(answer{2});      
         else
            load([short_fname,'_timesave.mat']);
        end
    end
    timestamp_fname = [analyze_times_directory,short_fname,'_timesave'];
   %% ali turned this off %% save(timestamp_fname,'start_time','end_time');
    %Log start and end times
    time_frames_temp = horzcat(file_to_analyze,num2cell(start_time),num2cell(end_time));
    if I == 1
        time_frames = time_frames_temp;
    else
        time_frames = vertcat(time_frames,time_frames_temp);
    end
    disp([file_to_analyze,' start and end times have been Logged.']);
    
    clear End_Temp_Array
    clear End_Temp 
end


%% % ------------------------- Load Files ------------------------- % %%

% File Location
cd(data_dir)    

% Preallocate Waveform_Index cell array
flds = [{'TS'}, {'Dat'}, {'SmpRt'}, {'B2V'}, {'File'}];
Waveform_Index = cell(file_total,length(flds));
       
for i=1:file_total
    file_to_analyze = cellstr([data_dir,files(i).name]);
    Waveform_Index(i,5) = file_to_analyze;
    % Get time stamps, sample points (i.e., time-series amplitude) and file
    % header
    [Waveform_Index{i,1}, Waveform_Index{i,2}, headr] = ...
        Nlx2MatCSC(file_to_analyze{1}, [1 0 0 0 1], 1, 1, []);

    % Parse header data
    Waveform_Index{i,3} = str2double(char(regexp(headr{14}, '(?<=\s)\S*', 'match')));
    Waveform_Index{i,4} = str2double(char(regexp(headr{16}, '(?<=\s)\S*', 'match')));

    % Convert data to 1D array
    Waveform_Index{i,2} = reshape(Waveform_Index{i,2},1,[]);
    % Convert from bit to voltage
    Waveform_Index{i,2} = Waveform_Index{i,2} * Waveform_Index{i,4};

    % Fill in ts values
    ts = repmat(Waveform_Index{i,1}, 512, 1) + ...
        repmat((linspace(0,511/Waveform_Index{i,3},512).*10^6)', 1, length(Waveform_Index{i,1}));
    Waveform_Index{i,1} = reshape(ts,1,[]);

    % Get Waveform of File
    W_Temp = resample(Waveform_Index{i,2},1000,Waveform_Index{i,3});
    Waveform(:,i) = W_Temp';

end
predata = cell2struct(Waveform_Index', flds, 1);
    
clearvars W_Temp
    
for iFile = 1:file_total %should loop through timestamped session folders   
    %% % ------------------------- File Analysis ------------------------- % %%
    
    %convert time in seconds to samples
    %ajw 5-9-13
    samples = [time_frames{iFile,2}*srate:time_frames{iFile,3}*srate];
    for channel=1
        %% % ------------------------- Colin Edit ------------------------ % %%
        
        % Colin 3/31/14 Modify code to set power threshold based on full
        % input file
        colins_edit = 1;
        if colins_edit
            % note: This section calculates power threshold but does not
            % account for clipping that might drive Pepisode higher
            [phase,pow]=multiphasevec2(freqs,Waveform(:,iFile)',srate,wavelet_cycles);
            zero_idx = find(pow ==0);  % find 0 values in pow
            pow_nan = pow;
            pow_nan(zero_idx) = NaN;  % convert 0 values in pow to NaN
            %%% from eeg toolbox calcPepisode
            %Blog = log10(double(B));
            Blog = log10(double(pow_nan));  % ajw modified
            Blog(zero_idx) = 0;
            Pm = mean(Blog,2);
            [all,R2] = chi_squarefit(freqs,Pm);
            all = all';
            %all = chi_squarefit_old(freqs,Pm)';
            % set the threshold
            pthresh = all(:,951);
            clearvars Blog all Pm pow_nan zero_idx
            % pthresh is now calculated, next lets find an amplitude threshold
            % to replace clipped values with NaNs
            [~,Athresh] = hist(Waveform(:,iFile)',100);
            % Athresh is going to serve as basis for removing clipping.
            nanclip = 0;
            if nanclip    
            else
                Athresh(1) = -Inf;
                Athresh(100) = Inf;
            end
        end
        
        
        %% % ------------------------- Singal Analysis ------------------------ % %%
        
        %Signal
        Signal = Waveform(samples,iFile);
        %         figure
        %         plot(Signal)
        %         tt = linspace(1,10,length(Signal));
        %         Signal = sin(2.*7.7.*pi.*tt)';
        %         hold on
        %         plot(Signal,'r')
        
        %P_Episode
        if colins_edit
            [~,~,~,~,Binary_matrix] = CK_Pepisode(Signal',freqs,wavelet_cycles,duration,srate,pthresh,Athresh);%,power_thresh)
        else
            [pow,phase,P_episode_pow,P_episode_phase,Binary_matrix] = AW_Pepisode(Signal',freqs,wavelet_cycles,duration,srate);%,power_thresh)
        end
        %
        
        %percentage time
        Percenttime(:,channel,iFile) = (sum(Binary_matrix,2)./size(Binary_matrix,2)).*100; %2 is the sum across the second dimension
        
        %added ajw 6-12-13
        tmp_freq_bin_pct_time = sum(Binary_matrix,1);
        tmp_freq_bin_pct_time(tmp_freq_bin_pct_time>1) = 1;
        freq_bin_pct_time(iFile,1) = (sum(tmp_freq_bin_pct_time)./length(tmp_freq_bin_pct_time))*100;
        clear tmp*
    end
    clear samples
end