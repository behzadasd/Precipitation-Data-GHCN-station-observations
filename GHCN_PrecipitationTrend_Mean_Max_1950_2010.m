%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mean and Maximum Precipitation Trends  %%%
%%%        GHCN Precipitation Data         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh, Ph.D.           %%%
%%% University of Pennsylvania      %%%
%%% basadieh@sas.upenn.edu          %%%
%%% github.com/behzadasd            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
clear;
clc;

%feature('UseGenericOpengl', 1); %%% Sets OpenGL to GenericOpengl to prevent the texts in plots being upside-down (due to a bug in OpenGL)
load Global_Mean_Temp_GISS_1901_2010
Global_Mean_Temp_GISS_1950_2010=Global_Mean_Temp_GISS_1901_2010(50:end,:); % for 1951 to 2010
Global_Mean_Temp_GISS_1950_2010(:,1)=[]; % Eliminates the 1st column which is the number of year for the corresponding temperature data

Continent_Names = {'Global';'North America'; 'South America'; 'Europe'; 'Oceania'; 'Africa'; 'Asia'; 'India'};

dir=[pwd '\']; % Current Directory Path
dir_data_in=[dir 'Matlab Data 20years\']; % Directory to raed raw data from
dir_out_fig=[pwd '\Figures\']; % Directory to save Figures and Maps
dir_mat_out=[pwd '\Results\'];

Days_M=[31; 28; 31; 30; 31; 30; 31; 31; 30; 31; 30;  31];
Year_first=1950;
Year_last=2010;
yrs_no=Year_last-Year_first+1; % Number of the years involved in the calculations
%yrs_no=61; % Number of the years involved in the calculations

min_No_yrs_prc=0.8; % mimimum percentage of days of the year with data to proceed with that year
min_NO_st_d=30; % Minimum number of available data for the station to have a reliable calculation

Input_File_Names = ls( fullfile(dir_data_in,'*.mat') ); % Lists the .mat files names in the specified directory
n_stations=size(Input_File_Names,1);

All_Lat=NaN(n_stations,1);
All_Lon=NaN(n_stations,1);
All_Elv=NaN(n_stations,1);
All_Station_ID=repmat(repmat(' ',1,11), n_stations,1); %creates empty chars with 11 empty spaces for the rows equal as the number of stations
All_Station_Names=repmat(repmat(' ',1,30), n_stations,1); %creates empty chars with 11 empty spaces for the rows equal as the number of stations

All_PRCP_annual=NaN(n_stations,yrs_no);
All_MaxPRCP_annual=NaN(n_stations,yrs_no);

for st=1:n_stations
    
    dir_data_name=[dir_data_in Input_File_Names(st,:)];  % Directory and Name of the Data to be read
    load (dir_data_name) % Loading the raw GCM Data
    
    All_Station_ID(st,:)=Station_ID;
    All_Station_Names(st,:)=Station_Name;
    All_Lat(st,1)=Lat;
    All_Lon(st,1)=Lon;
    All_Elv(st,1)=Elv;
    
    PRCP_daily=PRCP_daily * 0.1; % GHCN-Daily Precipittion unit is [tenths of mm]
    PRCP_daily(PRCP_daily<0)=NaN; % To eliminate missing data that are recorded as large negeative numbers
    
    %%% Creating the Precipitation time series vector for the station %%%
    PRCP_Long_Time_Series=NaN( yrs_no *365,1); % 70-yr precipitation time series of the station (each year is assumed 365-days fixed - No leap year consideration )
    %Date_Long_Time_Series=NaN( (No_decs*10) *365,1);
    for i=1:size(PRCP_daily,1)
        
        if ( str2double( Time(i,1:4) ) >= Year_first ) && ( str2double( Time(i,1:4) ) <= Year_last ) % checks to make sure the read month is in desired time interval
            
            dys_n= Days_M(str2double( Time(i,5:6) )); % Number of the days in the corresponding month to be read
            mnth_no=str2double( Time(i,5:6) ); % Count of the month
            yr_no=(str2double( Time(i,1:4) ) - Year_first ) +1; % Count of the year
            
            month_prcp=PRCP_daily(i,1:dys_n)';
            
            array_strt = (yr_no-1) *365 + sum(Days_M(1:mnth_no-1,1) ) +1 ; % The element number from which the monthly data should be fit in the long time series
            
            PRCP_Long_Time_Series( array_strt:array_strt+dys_n-1)= month_prcp;
            %Date_Long_Time_Series ( array_strt:array_strt+dys_n-1) = str2double( Time(i,1:6) ) ;
            
        end
        
    end
    
    %%% Calculating the Anuual-Mean daily precipitation value for each year for the station %%%
    for yr=1:yrs_no
        
        PRCP_year=PRCP_Long_Time_Series( (yr-1)*365+1 : yr*365 ,1); % Precipitation time series for the selected decade
        
        NaN_Check=~isnan(PRCP_year);
        if sum(NaN_Check,1) >= 365 * min_No_yrs_prc % That means the station has the mimimum percentage of daily precipitation data for this year
            
            All_PRCP_annual(st,yr)=nanmean(PRCP_year);
            All_MaxPRCP_annual(st,yr)=nanmax(PRCP_year);
        end
        
    end
    
    if sum(All_PRCP_annual(st,:)==0) >= 1 % That means this station has a year with PRCP=0, so this station is not suitable for RWHS
        All_PRCP_annual(st,:)=NaN;
        All_MaxPRCP_annual(st,:)=NaN;
    end
    
    text=['GHCN Daily Precipitation - St. Name= ', Station_ID, ' *** File ', num2str(st), ' / ', num2str(n_stations) ' Processed'];
    disp(text) % this part prints the current number of itteraion and its best solution and fitness at each iteration
        
end

[All_b_slope_PRCP_ann, All_b_prpc_temp_ann, All_PRCP_ave_ann, ~, ~, NO_st_d_PRCP_ann]=func_LinReg2_Trend_station(All_PRCP_annual, Global_Mean_Temp_GISS_1950_2010, min_NO_st_d);
[All_b_slope_MaxPRCP_ann, All_b_maxprpc_temp_ann, All_MaxPRCP_ave_ann, ~, ~, NO_st_d_MaxPRCP_ann]=func_LinReg2_Trend_station(All_MaxPRCP_annual, Global_Mean_Temp_GISS_1950_2010, min_NO_st_d);
All_b_prpc_temp_ann = All_b_prpc_temp_ann * 100;
All_b_maxprpc_temp_ann = All_b_maxprpc_temp_ann * 100;

%%% Fractional change in Precipitation using Natural Logarithm (Ln) %%%
[All_b_slope_LogPRCP_ann, All_LogPRCP_ave_ann, ~, ~, ~]=func_LinReg1_Trend_station(log(All_PRCP_annual), min_NO_st_d);
[All_b_slope_LogMaxPRCP_ann, All_LogMaxPRCP_ave_ann, ~, ~, ~]=func_LinReg1_Trend_station(log(All_MaxPRCP_annual), min_NO_st_d);
All_Frac_prcp_b_ann_1= ( exp (All_b_slope_LogPRCP_ann) -1 ) *100 ; % Relative Change in Annual-Average Mean Precipitation unisng the natural logarithm
All_Frac_maxprcp_b_ann_1= ( exp (All_b_slope_LogMaxPRCP_ann) -1 ) *100 ; % Relative Change in Annual-Average Maximum Precipitation unisng the natural logarithm

%%% Eliminating stations with NaN values %%%
for i=size(All_b_slope_PRCP_ann,1):-1:1
    if isnan(All_b_slope_PRCP_ann(i,1))
        
        All_PRCP_annual(i,:)=[];
        All_MaxPRCP_annual(i,:)=[];
        All_b_slope_PRCP_ann(i,:)=[];
        All_b_slope_MaxPRCP_ann(i,:)=[];
        All_b_prpc_temp_ann(i,:)=[];
        All_b_maxprpc_temp_ann(i,:)=[];
        All_PRCP_ave_ann(i,:)=[];
        All_MaxPRCP_ave_ann(i,:)=[];
        All_Frac_prcp_b_ann_1(i,:)=[];
        All_Frac_maxprcp_b_ann_1(i,:)=[];
        All_Station_ID(i,:)=[];
        All_Station_Names(i,:)=[];
        All_Lat(i,:)=[];
        All_Lon(i,:)=[];
        All_Elv(i,:)=[];
        NO_st_d_PRCP_ann(i,:)=[];
        NO_st_d_MaxPRCP_ann(i,:)=[];
        All_b_slope_LogMaxPRCP_ann(i,:)=[];
        All_b_slope_LogPRCP_ann(i,:)=[];
        All_LogMaxPRCP_ave_ann(i,:)=[];
        All_LogPRCP_ave_ann(i,:)=[];
        
    end
end

%%% Fractional change in Precipitation using b/<P> %%%
All_Frac_prcp_b_ann_2=NaN(size(All_b_slope_PRCP_ann,1),1); %  Fractional change in Precipitation b/<P> for each station
for i=1:size(All_b_slope_PRCP_ann,1)
    if ~isnan(All_b_slope_PRCP_ann(i,1)) && ~isnan(All_PRCP_ave_ann(i,1))
        if All_PRCP_ave_ann(i,1)==0 % to ignore denominator of the fraction being zero leading to variable being NaN
            All_Frac_prcp_b_ann_2(i,1)=0;
        else
            All_Frac_prcp_b_ann_2(i,1)=(All_b_slope_PRCP_ann(i,1)/All_PRCP_ave_ann(i,1)) * 100;
        end
    end
end

All_Frac_maxprcp_b_ann_2=NaN(size(All_b_slope_MaxPRCP_ann,1),1); %  Fractional change in Precipitation b/<P> for each station
for i=1:size(All_b_slope_MaxPRCP_ann,1)
    if ~isnan(All_b_slope_MaxPRCP_ann(i,1)) && ~isnan(All_MaxPRCP_ave_ann(i,1))
        if All_MaxPRCP_ave_ann(i,1)==0 % to ignore denominator of the fraction being zero leading to variable being NaN
            All_Frac_maxprcp_b_ann_2(i,1)=0;
        else
            All_Frac_maxprcp_b_ann_2(i,1)=(All_b_slope_MaxPRCP_ann(i,1)/All_MaxPRCP_ave_ann(i,1)) * 100;
        end
    end
end

%%% Precipitation-Weighted Global Average %%%

All_Results_prcp_weighted(1,1) = nansum(nansum( All_b_prpc_temp_ann .* All_PRCP_ave_ann )) / nansum(nansum( All_PRCP_ave_ann )) ;
All_Results_prcp_weighted(1,2) = nansum(nansum( All_b_maxprpc_temp_ann .* All_MaxPRCP_ave_ann )) / nansum(nansum( All_MaxPRCP_ave_ann )) ;

%%% Trend Averaging by Continent %%%

[C_Ave_PRCP_annual, C_Ave_MaxPRCP_annual]=func_TrendAve_2var_station(All_PRCP_ave_ann, All_MaxPRCP_ave_ann, All_Lat, All_Lon, NO_st_d_PRCP_ann, min_NO_st_d);
[C_Ave_b_slope_PRCP_ann, C_Ave_b_slope_MaxPRCP_ann]=func_TrendAve_2var_station(All_b_slope_PRCP_ann, All_b_slope_MaxPRCP_ann, All_Lat, All_Lon, NO_st_d_PRCP_ann, min_NO_st_d);
[C_Ave_Frac_prcp_b_ann_2, C_Ave_Frac_maxprcp_b_ann_2]=func_TrendAve_2var_station(All_Frac_prcp_b_ann_2, All_Frac_maxprcp_b_ann_2, All_Lat, All_Lon, NO_st_d_PRCP_ann, min_NO_st_d);
[C_Ave_Frac_prcp_b_ann_1, C_Ave_Frac_maxprcp_b_ann_1]=func_TrendAve_2var_station(All_Frac_prcp_b_ann_1, All_Frac_maxprcp_b_ann_1, All_Lat, All_Lon, NO_st_d_PRCP_ann, min_NO_st_d);
[C_Ave_b_prpc_temp_ann, C_Ave_b_maxprpc_temp_ann]=func_TrendAve_2var_station(All_b_prpc_temp_ann, All_b_maxprpc_temp_ann, All_Lat, All_Lon, NO_st_d_PRCP_ann, min_NO_st_d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
All_Results=NaN(8,10); % Stores the global average of results of " PRCP (mm/day), b (mm/day/yr), frac_b (% / yr)(b/<P>), frac_b (% / yr)(Ln), b_per_t (% / K) " for Mean and Max Precipitation (overal 10 variable)
All_Results(:,1)=C_Ave_PRCP_annual;
All_Results(:,2)=C_Ave_MaxPRCP_annual;
All_Results(:,3)=C_Ave_b_slope_PRCP_ann;
All_Results(:,4)=C_Ave_b_slope_MaxPRCP_ann;
All_Results(:,5)=C_Ave_Frac_prcp_b_ann_2;
All_Results(:,6)=C_Ave_Frac_maxprcp_b_ann_2;
All_Results(:,7)=C_Ave_Frac_prcp_b_ann_1;
All_Results(:,8)=C_Ave_Frac_maxprcp_b_ann_1;
All_Results(:,9)=C_Ave_b_prpc_temp_ann;
All_Results(:,10)=C_Ave_b_maxprpc_temp_ann;

excel_out_name=[pwd '\GHCN-Daily Max-Mean Precipitation 1950-2010 - Results Statistics summary.xls'];
variables1={'Ave PRCP [mm/day]', 'Ave PRCP [mm]', 'b [mm/day/yr]', 'b [mm/yr]', 'frac_b [% / yr] (b/<P>)', 'frac_b [% / yr]  (b/<P>)', 'frac_b [% / yr] (ln(P))', 'frac_b [% / yr]  (ln(P))', 'b_per_t [% / K]','b_per_t [% / K]'};
variables2={'Mean PRCP', 'Max PRCP', 'Mean PRCP', 'Max PRCP', 'Mean PRCP', 'Max PRCP', 'Mean PRCP', 'Max PRCP', 'Mean PRCP', 'Max PRCP'};
xlswrite(excel_out_name, variables1, 'Sheet1', 'B3')
xlswrite(excel_out_name, variables2, 'Sheet1', 'B4')
xlswrite(excel_out_name, cellstr('GHCN-Daily 1950-2010'), 'Sheet1', 'A1')
xlswrite(excel_out_name, cellstr(Continent_Names), 'Sheet1', 'A5')
xlswrite(excel_out_name, All_Results, 'Sheet1', 'B5')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Percentile=prctile(All_PRCP_ave_ann, [10 20 30 40 50 60 70 80 90 100]);
Percentiled_PRCP_ave_ann=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);
Percentiled_MaxPRCP_ave_ann=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);
Percentiled_b_slope_PRCP_ann=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);
Percentiled_b_slope_MaxPRCP_ann=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);
Percentiled_b_prpc_temp_ann=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);
Percentiled_b_maxprpc_temp_ann=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);
Percentiled_Frac_prcp_b_ann_1=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);
Percentiled_Frac_maxprcp_b_ann_1=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);
Percentiled_Frac_prcp_b_ann_2=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);
Percentiled_Frac_maxprcp_b_ann_2=NaN( ceil( (size(All_PRCP_ave_ann,1)/10)+10 ) ,10);

ct=zeros(1,10);

for i=size(All_PRCP_ave_ann,1):-1:1
    
    if All_PRCP_ave_ann(i,1) > 0 && All_PRCP_ave_ann(i,1) <= Percentile(1,1)
        
        ct(1,1)=ct(1,1)+1;
        Percentiled_PRCP_ave_ann(ct(1,1),1)=All_PRCP_ave_ann(i,1);
        Percentiled_MaxPRCP_ave_ann(ct(1,1),1)=All_MaxPRCP_ave_ann(i,1);
        Percentiled_b_slope_PRCP_ann(ct(1,1),1)=All_b_slope_PRCP_ann(i,1);
        Percentiled_b_slope_MaxPRCP_ann(ct(1,1),1)=All_b_slope_MaxPRCP_ann(i,1);
        Percentiled_b_prpc_temp_ann(ct(1,1),1)=All_b_prpc_temp_ann(i,1);
        Percentiled_b_maxprpc_temp_ann(ct(1,1),1)=All_b_maxprpc_temp_ann(i,1);
        Percentiled_Frac_prcp_b_ann_1(ct(1,1),1)=All_Frac_prcp_b_ann_1(i,1);
        Percentiled_Frac_maxprcp_b_ann_1(ct(1,1),1)=All_Frac_maxprcp_b_ann_1(i,1);
        Percentiled_Frac_prcp_b_ann_2(ct(1,1),1)=All_Frac_prcp_b_ann_2(i,1);
        Percentiled_Frac_maxprcp_b_ann_2(ct(1,1),1)=All_Frac_maxprcp_b_ann_2(i,1);
        
    end
    
    flag=0;
    for j=2:10
        
        if All_PRCP_ave_ann(i,1) > Percentile(1,j-1) && All_PRCP_ave_ann(i,1) <= Percentile(1,j)
            
            ct(1,j)=ct(1,j)+1;            
            Percentiled_PRCP_ave_ann(ct(1,j),j)=All_PRCP_ave_ann(i,1);
            Percentiled_MaxPRCP_ave_ann(ct(1,j),j)=All_MaxPRCP_ave_ann(i,1);
            Percentiled_b_slope_PRCP_ann(ct(1,j),j)=All_b_slope_PRCP_ann(i,1);
            Percentiled_b_slope_MaxPRCP_ann(ct(1,j),j)=All_b_slope_MaxPRCP_ann(i,1);
            Percentiled_b_prpc_temp_ann(ct(1,j),j)=All_b_prpc_temp_ann(i,1);
            Percentiled_b_maxprpc_temp_ann(ct(1,j),j)=All_b_maxprpc_temp_ann(i,1);
            Percentiled_Frac_prcp_b_ann_1(ct(1,j),j)=All_Frac_prcp_b_ann_1(i,1);
            Percentiled_Frac_maxprcp_b_ann_1(ct(1,j),j)=All_Frac_maxprcp_b_ann_1(i,1);
            Percentiled_Frac_prcp_b_ann_2(ct(1,j),j)=All_Frac_prcp_b_ann_2(i,1);
            Percentiled_Frac_maxprcp_b_ann_2(ct(1,j),j)=All_Frac_maxprcp_b_ann_2(i,1);
            flag=1;
            
        end
        if flag ==1
            break
        end
        
    end
    
end

Percentiled_Ave_PRCP_ave_ann=nanmean(Percentiled_PRCP_ave_ann);
Percentiled_Ave_MaxPRCP_ave_ann=nanmean(Percentiled_MaxPRCP_ave_ann);
Percentiled_Ave_b_slope_PRCP_ann=nanmean(Percentiled_b_slope_PRCP_ann);
Percentiled_Ave_b_slope_MaxPRCP_ann=nanmean(Percentiled_b_slope_MaxPRCP_ann);
Percentiled_Ave_b_prpc_temp_ann=nanmean(Percentiled_b_prpc_temp_ann);
Percentiled_Ave_b_maxprpc_temp_ann=nanmean(Percentiled_b_maxprpc_temp_ann);
Percentiled_Ave_Frac_prcp_b_ann_1=nanmean(Percentiled_Frac_prcp_b_ann_1);
Percentiled_Ave_Frac_maxprcp_b_ann_1=nanmean(Percentiled_Frac_maxprcp_b_ann_1);
Percentiled_Ave_Frac_prcp_b_ann_2=nanmean(Percentiled_Frac_prcp_b_ann_2);
Percentiled_Ave_Frac_maxprcp_b_ann_2=nanmean(Percentiled_Frac_maxprcp_b_ann_2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('All_Results') % Saves all the variables in the workspace in the current directory

GHCN_C_Ave_PRCP_Mean=C_Ave_PRCP_annual';
GHCN_C_Ave_PRCP_Max=C_Ave_MaxPRCP_annual';
GHCN_C_Ave_b_slope_PRCP_Mean=C_Ave_b_slope_PRCP_ann';
GHCN_C_Ave_b_slope_PRCP_Max=C_Ave_b_slope_MaxPRCP_ann';
GHCN_C_Ave_Frac_prcp_b_2_Mean=C_Ave_Frac_prcp_b_ann_2';
GHCN_C_Ave_Frac_prcp_b_2_Max=C_Ave_Frac_maxprcp_b_ann_2';
GHCN_C_Ave_Frac_prcp_b_1_Mean=C_Ave_Frac_prcp_b_ann_1';
GHCN_C_Ave_Frac_prcp_b_1_Max=C_Ave_Frac_maxprcp_b_ann_1';
GHCN_C_Ave_b_prcp_temp_Mean=C_Ave_b_prpc_temp_ann';
GHCN_C_Ave_b_prcp_temp_Max=C_Ave_b_maxprpc_temp_ann';
  

%%%%%%%%%%%%%%%%%%%
%%% Global Maps %%%
dir_out_name_fig_tag='_GHCN_1950-2010';
x_limit=[-180 180]; y_limit=[-90 90];

dir_out_name_fig=[dir_out_fig 'PRCP_Ave_' dir_out_name_fig_tag];
Variable=All_PRCP_ave_ann;
limits = max(abs(quantile(Variable(:), [0.05 0.95])));
var_limit=[0 limits];
title_text={'Annual-Averaged Daily Precipitation [mm.day^-^1]'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'};
[~]=func_GeoshowMapScatter_station(Variable, All_Lat, All_Lon, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

dir_out_name_fig=[dir_out_fig 'MaxPRCP_Ave_' dir_out_name_fig_tag];
Variable=All_MaxPRCP_ave_ann;
limits = max(abs(quantile(Variable(:), [0.05 0.95])));
var_limit=[0 limits];
title_text={'Annual-Maximum Daily Precipitation [mm]'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'};
[~]=func_GeoshowMapScatter_station(Variable, All_Lat, All_Lon, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

dir_out_name_fig=[dir_out_fig 'b_PRCP_' dir_out_name_fig_tag];
Variable=All_b_slope_PRCP_ann;
limits = max(abs(quantile(Variable(:), [0.05 0.95])));
var_limit=[0 limits];
title_text={'Slope of Change (b) in Annual-Averaged Daily Precipitation [mm.day^-^1/year]'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'};
[~]=func_GeoshowMapScatter_station(Variable, All_Lat, All_Lon, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

dir_out_name_fig=[dir_out_fig 'b_MaxPRCP_' dir_out_name_fig_tag];
Variable=All_b_slope_MaxPRCP_ann;
limits = max(abs(quantile(Variable(:), [0.05 0.95])));
var_limit=[0 limits];
title_text={'Slope of Change (b) in Annual-Maximum Daily Precipitation [mm/year]'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'};
[~]=func_GeoshowMapScatter_station(Variable, All_Lat, All_Lon, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

dir_out_name_fig=[dir_out_fig 'b_PRCP_Frac_' dir_out_name_fig_tag];
Variable=All_Frac_prcp_b_ann_1;
%limits = max(abs(quantile(Variable(:), [0.05 0.95])));
%var_limit=[-limits limits];
var_limit=[-1 1];
title_text={'Relative Change (b) in Annual-Averaged Daily Precipitation [%.year^-^1]'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'};
[~]=func_GeoshowMapScatter_station(Variable, All_Lat, All_Lon, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

dir_out_name_fig=[dir_out_fig 'b_MaxPRCP_Frac_' dir_out_name_fig_tag];
Variable=All_Frac_maxprcp_b_ann_1;
%limits = max(abs(quantile(Variable(:), [0.05 0.95])));
%var_limit=[-limits limits];
var_limit=[-1 1];
title_text={'Relative Change (b) in Annual-Maximum Daily Precipitation [%.year^-^1]'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'};
[~]=func_GeoshowMapScatter_station(Variable, All_Lat, All_Lon, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

dir_out_name_fig=[dir_out_fig 'b_prcp_temp_' dir_out_name_fig_tag];
Variable=All_b_prpc_temp_ann;
%limits = max(abs(quantile(Variable(:), [0.05 0.95])));
%var_limit=[-limits limits];
var_limit=[-80 80];
title_text={'Relative Change (b) in Annual-Averaged Daily Precipitation per K Global Warming [%.K^-^1]'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'};
[~]=func_GeoshowMapScatter_station(Variable, All_Lat, All_Lon, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

dir_out_name_fig=[dir_out_fig 'b_maxprcp_temp_' dir_out_name_fig_tag];
Variable=All_b_maxprpc_temp_ann;
%limits = max(abs(quantile(Variable(:), [0.05 0.95])));
%var_limit=[-limits limits];
var_limit=[-80 80];
title_text={'Relative Change (b) in Annual-Maximum Daily Precipitation per K Global Warming [%.K^-^1]'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'};
[~]=func_GeoshowMapScatter_station(Variable, All_Lat, All_Lon, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

%%%%%%%%%%%%%%%%%%%%%
%%% Scatter Plots %%%
dir_out_name_fig=[dir_out_fig 'Scatter_b_P_MeanMax_relative_' dir_out_name_fig_tag];
h1=scatter(All_Frac_prcp_b_ann_1, All_Frac_maxprcp_b_ann_1, 'o','MarkerFaceColor',[0 0.6 0.6],'MarkerEdgeColor',[0.02 0.02 0.02],'LineWidth',0.5, 'SizeData',20);
title({'Sactter Plot of Relative Change in Mean Precipitation vs. Max Precipitation - Real Time Precipitation'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'});
xlabel({'Relative Change in Mean Daily Precipitation'; '[%.decade^-^1]'});
ylabel({'Relative Change in Maximum Daily Precipitation'; '[%.decade^-^1]'});
set(gca,'FontSize',26) % Axis Numbers and rages Font
set(findall(gcf,'type','text'),'FontSize',22,'FontWeight', 'normal') % Text (title and axis-lable) font
axis square
lsline % Adds least-square line
%refline(1,-5)
xlim([-20 20]); ylim([-20 20]);
hold on
plot([-50,100 * cosd(10*4.5)], [-50,100 * sind(10*4.5)],'LineWidth',1,'Color','black'); % identity (1by1) line
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
saveas(h1,dir_out_name_fig)
saveas(h1,dir_out_name_fig,'png')


%%%%%%%%%%%%%
%%% Plots %%%
Plot_xtick={'0-10th', '10-20th', '20-30th', '30-40th', '40-50th', '50-60th', '60-70th', '70-80th', '80-90th', '90-100th'};

dir_out_name_fig=[dir_out_fig 'Plot_b_PRCP_Frac_' dir_out_name_fig_tag];
h1=plot((1:10), Percentiled_Ave_Frac_prcp_b_ann_2, 'b-d', (1:10), Percentiled_Ave_Frac_maxprcp_b_ann_2, 'r-o', 'LineWidth',4, 'MarkerSize',20);
ylim([-0.02 0.18]);
%h1=plot((1:10), Percentiled_Ave_Frac_prcp_b_ann_2, 'b-d','LineWidth',4, 'MarkerSize',20, 'MarkerEdgeColor',[0 0.6 0.6], 'MarkerFaceColor',[0 0.6 0.6]);
set(gca,'XTickLabel',Plot_xtick)
legend('Mean Precipitation', 'Maximum Precipitation')
title({'Relative Change (b) in Mean and Maximum Precipitation [%.year^-^1] by precipitation percentile'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'})
xlabel({'Annual Precipitation Percentiles'});
ylabel({'Relative Change (b) in Precipitation'; '[%.year^-^1]'});
set(gca,'FontSize',26) % Axis Numbers and rages Font
set(findall(gcf,'type','text'),'FontSize',22,'FontWeight', 'normal') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
saveas(h1,dir_out_name_fig)
saveas(h1,dir_out_name_fig,'png')

dir_out_name_fig=[dir_out_fig 'Plot_b_prcp_temp_' dir_out_name_fig_tag];
h1=plot((1:10), Percentiled_Ave_b_prpc_temp_ann, 'b-d', (1:10), Percentiled_Ave_b_maxprpc_temp_ann, 'r-o', 'LineWidth',4, 'MarkerSize',20);
ylim([0 20]);
%h1=plot((1:10), Percentiled_Ave_Frac_prcp_b_ann_2, 'b-d','LineWidth',4, 'MarkerSize',20, 'MarkerEdgeColor',[0 0.6 0.6], 'MarkerFaceColor',[0 0.6 0.6]);
set(gca,'XTickLabel',Plot_xtick)
legend('Mean Precipitation', 'Maximum Precipitation')
title({'Relative Change (b) in Mean and Maximum Precipitation per K Global Warming [%.K^-^1] by precipitation percentile'; 'GHCN-Daily Observational Data - Max/Mean Precipitation 1950-2010'})
xlabel({'Annual Precipitation Percentiles'});
ylabel({'Relative Change (b) in Precipitation per K Global Warming'; '[%.K^-^1]'});
set(gca,'FontSize',26) % Axis Numbers and rages Font
set(findall(gcf,'type','text'),'FontSize',22,'FontWeight', 'normal') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% saveas(h1,dir_out_name_fig)
% saveas(h1,dir_out_name_fig,'png')


toc;

