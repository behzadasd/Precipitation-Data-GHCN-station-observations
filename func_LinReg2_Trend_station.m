function [b_slope, b_data1_data2, Data_ave, All_stast, Reg_Significance, NO_st_d]=func_LinReg2_Trend_station(Data1, Data2, min_NO_st_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Linear Regression Trend Analysis         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data1 = PRCP
% Data2 = Temp
% b_data1_data2= b_prcp_temp = Relative Change in Precipitation Regarding the Global Mean Temperature
stations_n=size(Data1,1); % Total number of the stations
yrs_n=size(Data1,2); % Number of the Years that data is available for

Strt_yr=zeros(stations_n,1); % The count of year that the data is available from, for each station
NO_st_d=zeros(stations_n,1); % Number of the data available for each station

b_slope=NaN(stations_n,1);
b_data1_data2=NaN(stations_n,1);
bint_05=NaN(stations_n,2,2); % Intervals of b at 95% confidence (5% Significance Level)
bint_01=NaN(stations_n,2,2); % Intervals of b at 99% confidence (1% Significance Level)
All_stast=NaN(stations_n,4);
Reg_Significance=NaN(stations_n,2); % Linear Regression trend Significance at 1% and 5% level (1st and 2nd layer, respectively)
All_Z1_05=NaN(stations_n,1); % Z1= b / delta(b) at 5% Level
All_Z1_01=NaN(stations_n,1); % Z1= b / delta(b) at 1% Level
Data_ave=NaN(stations_n,1); % Average Data of each station

for St=1:stations_n
        
        for i=1:yrs_n
            if ~isnan(Data1(St,i))
                Strt_yr(St,1)=i; % Gives the count of year that the data of this grid has been started from
                break
            end
        end
        
        if Strt_yr(St,1)>0 % That means this grid has at least 1 data available
            counter=0;
            for i=Strt_yr(St,1):yrs_n
                if ~isnan(Data1(St,i))
                    counter=counter+1;
                end
            end
            NO_st_d(St,1)=counter; % Number of the data available for this grid
        end
        
        if NO_st_d(St,1)>= 1 % That means this station has at least 1 data available
            Reg_Significance(St,:)=0;
        end
        
        if NO_st_d(St,1)>= min_NO_st_d %%% Minimum number of available data for the station to have a reliable calculation %%%
            
            t_vec=zeros(NO_st_d(St,1),1);
            counter_a=1;
            counter_b=1;
            for i=1:yrs_n
                if ~isnan(Data1(St,i))
                    t_vec(counter_b,1)=counter_a;
                    counter_b=counter_b+1;
                end
                counter_a=counter_a+1;
            end
            
            help_Data=zeros(NO_st_d(St,1),1); % Gathers all non-NaN values of the current station in a vector
            help_Data_temp=zeros(NO_st_d(St,1),1); % Creates a Temp vector for the years with PRCP data available
            counter=1;
            for i=1:yrs_n
                if ~isnan(Data1(St,i))
                    help_Data(counter,1)=Data1(St,i);
                    help_Data_temp(counter,1)=Data2(i,1);
                    counter=counter+1;
                end
            end
            
            x_vec=[ones(NO_st_d(St,1),1) t_vec];
            temp_vec=[ones(NO_st_d(St,1),1) help_Data_temp];
            
            %%% Relative Change in Precipitation Regarding the Global Mean Temperature %%%
            
            bt=regress(log(help_Data), temp_vec);
            b_data1_data2(St,1)=exp(bt(2,1))-1;
                       
            %%% 5% Level Calculations %%%
            [b,bint,~,~,stats] = regress(help_Data,x_vec);
            
            b_slope(St,1)=b(2,1);
            bint_05(St,:,1)=bint(1,:);
            bint_05(St,:,2)=bint(2,:);
            All_stast(St,:)=stats;
            Data_ave(St,1)=mean(help_Data);
            
            All_Z1_05(St,1)=b(2,1)/((abs(bint(2,2)-bint(2,1)))/2);
            
            if ( stats(1,3) < 0.05 && b(2,1)>=0 )
                Reg_Significance(St,2)=1;
            elseif ( stats(1,3) < 0.05 && b(2,1)<0 )
                Reg_Significance(St,2)=-1;
            end
            
            %%% 1% Level Calculations %%%
            [b,bint,~,~,stats] = regress(help_Data,x_vec,0.01);
            
            bint_01(St,:,1)=bint(1,:);
            bint_01(St,:,2)=bint(2,:);
            
            All_Z1_01(St,1)=b(2,1)/((abs(bint(2,2)-bint(2,1)))/2);
            
            if ( stats(1,3) < 0.01 && b(2,1)>=0 )
                Reg_Significance(St,1)=1;
            elseif ( stats(1,3) < 0.01 && b(2,1)<0 )
                Reg_Significance(St,1)=-1;
            end
            
        end

end


end

