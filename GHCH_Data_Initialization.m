%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    GHCN Daily Precipitation Data Load   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh, Ph.D.           %%%
%%% University of Pennsylvania      %%%
%%% basadieh@sas.upenn.edu          %%%
%%% github.com/behzadasd            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
clear;
clc;

%feature('UseGenericOpengl', 1); %%% Sets OpenGL to GenericOpengl to prevent the texts in plots being upside-down (due to a bug in OpenGL)

variable_id = 'PRCP';
min_months= 20 * 12; % minimum number of months with available precipitation data

start_year = 1850;
end_year = 2012;


dir=[pwd '\']; % Current Directory Path
dir_data_in=[dir 'ghcnd_all\']; % Directory to read raw data from
dir_mat_out=[pwd '\Matlab Data 20years\'];

Input_File_Names = ls( fullfile(dir_data_in,'*.dly') ); % Lists the .dly files names in the specified directory

n_stations=size(Input_File_Names,1);

station_filename = [pwd '\ghcnd-stations.txt'];
fid = fopen(station_filename, 'r');
Station_data = fscanf(fid,'%86c \n', [86, Inf]); %text in fixed-width columns; size: 86*n_stations
fclose(fid);
Station_data=Station_data';   Station_data(1,:)=[];
%A=Station_data(1:100,:);

PRCP_GHCN_1951_2014=NaN(2, 23376);


for s = 1:n_stations
    
    %read the appropriate file
    filename = [dir_data_in   Input_File_Names(s,:)];
    
    fid = fopen(filename, 'r');
    All_data = fscanf(fid,'%269c \n', [269, Inf]); %21 header characters plus 31*8 data+flags
    fclose(fid);
    
    All_data=All_data';
    
    
    for i=size(All_data,1):-1:1
        if ~strncmp(All_data(i,18:21), variable_id, 4) % All_data(i,18:21) represents the ELEMENT or the name of the variable
            All_data(i,:)=[]; % If the name of the variable is not the same as the target variable (here PRCP), delets that line
        end
    end
    
    if size(All_data,1) > min_months % Proceed if the station has the minimum number of months with data
        
        Station_ID=All_data(1,1:11);
        PRCP_daily=NaN(size(All_data,1),31);
        Time = repmat('      ',size(PRCP_daily,1),1);
        
        for YM=1:size(PRCP_daily,1)
            for D=1:31
                
                PRCP_daily(YM,D)=str2double( All_data (YM , 22+(D-1)*8 : 22+(D-1)*8+4 ) );
                Time(YM,1:6)=All_data(YM,12:17);
            end
        end
        PRCP_daily(PRCP_daily==-9999)=NaN;
        
        
        for i=1:size(Station_data,1)
            
            if strncmp( Station_ID, Station_data(i,1:11), 11)
                st_no=i; % Station Number on the list (to read the lat, lon, elv and etc. for the station from that line)
                Lat=str2double( Station_data(st_no,13:20) );
                Lon=str2double( Station_data(st_no,22:30) );
                Elv=str2double( Station_data(st_no,32:37) );
                Station_Name=Station_data(st_no,42:71);
                
                break
            end
            
        end
        
        
        %Station_data(st_no,:)=[]; % Delets the station from the list after reading to make the search for the remaining stations faster
        
        file_mat_dir_out=[dir_mat_out  Station_ID '_GHCN_PRCP_Daily' '.mat']; % Directory and Name of .mat file to be saved after analysis
        
        save(file_mat_dir_out, 'PRCP_daily', 'Time', 'Lat', 'Lon', 'Elv', 'Station_Name', 'Station_ID', 'variable_id', 'st_no')
        disp(['Saved   ' Station_ID '   ***   File ' num2str(s) ' of total ' num2str(n_stations) ' Processed'])
        
    end
    
end



toc;











