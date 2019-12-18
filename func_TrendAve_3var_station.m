function [C_Ave_Var1, C_Ave_Var2, C_Ave_Var3]=func_TrendAve_3var_station(Var1, Var2, Var3, Lat, Lon, NO_st_d, min_NO_st_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Trend test results Averaging - Globally and by Continent   %%%

load Continents_Bounds;
% Columns are: Continents are North America, South America, Oceania, Europe1, Europe2, Asia 1-6, India1, India2, Africa 1, Africa 2
% India will be included in Asia average as well
% Rows are: Lat_Down, Lat_Up, Lon_Left, Lon_Right

stations_n=size(Var1,1); % Total number of the stations

sum_Var1=zeros(size(Continents_LatLong_Bound,1),1);
sum_Var2=zeros(size(Continents_LatLong_Bound,1),1);
sum_Var3=zeros(size(Continents_LatLong_Bound,1),1);

station_counter=zeros(size(Continents_LatLong_Bound,1),1); % number of stations for each continent range

for i=1:stations_n
    
    if NO_st_d(i,1)>= min_NO_st_d
        
        for j=1:size(sum_Var1,1)
            
            if (Lat(i,1) > Continents_LatLong_Bound(j,1)) && (Lat(i,1) <= Continents_LatLong_Bound(j,2))
                if (Lon(i,1) > Continents_LatLong_Bound(j,3)) && (Lon(i,1) <= Continents_LatLong_Bound(j,4))
                    
                    if ~isnan(Var1(i,1))
                        station_counter(j,1)=station_counter(j,1)+1;
                        sum_Var1(j,1)=sum_Var1(j,1)+Var1(i,1);
                    end
                    if ~isnan(Var2(i,1))
                        sum_Var2(j,1)=sum_Var2(j,1)+Var2(i,1);
                    end
                    if ~isnan(Var3(i,1))
                        sum_Var3(j,1)=sum_Var3(j,1)+Var3(i,1);
                    end
                    
                end
            end
            
        end
        
    end
    
end

% Global, North America, South America, Europe, Oceania, Africa, Asia, India
C_Ave_Var1=zeros(8,1);

C_Ave_Var1(1,1)=nanmean(Var1,1); % Global
C_Ave_Var1(2,1)=sum_Var1(1,1)/station_counter(1,1); % North America
C_Ave_Var1(3,1)=sum_Var1(2,1)/station_counter(2,1); % South America
C_Ave_Var1(4,1)=sum(sum_Var1(4:5,1))/sum(station_counter(4:5,1)); % Europe
C_Ave_Var1(5,1)=sum_Var1(3,1)/station_counter(3,1); % Oceania
C_Ave_Var1(6,1)=sum(sum_Var1(14:15,1))/sum(station_counter(14:15,1)); % Africa
C_Ave_Var1(7,1)=sum(sum_Var1(6:13,1))/sum(station_counter(6:13,1)); % Asia
C_Ave_Var1(8,1)=sum(sum_Var1(12:13,1))/sum(station_counter(12:13,1)); % India

C_Ave_Var2=zeros(8,1);

C_Ave_Var2(1,1)=nanmean(Var2,1); % Global
C_Ave_Var2(2,1)=sum_Var2(1,1)/station_counter(1,1); % North America
C_Ave_Var2(3,1)=sum_Var2(2,1)/station_counter(2,1); % South America
C_Ave_Var2(4,1)=sum(sum_Var2(4:5,1))/sum(station_counter(4:5,1)); % Europe
C_Ave_Var2(5,1)=sum_Var2(3,1)/station_counter(3,1); % Oceania
C_Ave_Var2(6,1)=sum(sum_Var2(14:15,1))/sum(station_counter(14:15,1)); % Africa
C_Ave_Var2(7,1)=sum(sum_Var2(6:13,1))/sum(station_counter(6:13,1)); % Asia
C_Ave_Var2(8,1)=sum(sum_Var2(12:13,1))/sum(station_counter(12:13,1)); % India

C_Ave_Var3=zeros(8,1);

C_Ave_Var3(1,1)=nanmean(Var3,1); % Global
C_Ave_Var3(2,1)=sum_Var3(1,1)/station_counter(1,1); % North America
C_Ave_Var3(3,1)=sum_Var3(2,1)/station_counter(2,1); % South America
C_Ave_Var3(4,1)=sum(sum_Var3(4:5,1))/sum(station_counter(4:5,1)); % Europe
C_Ave_Var3(5,1)=sum_Var3(3,1)/station_counter(3,1); % Oceania
C_Ave_Var3(6,1)=sum(sum_Var3(14:15,1))/sum(station_counter(14:15,1)); % Africa
C_Ave_Var3(7,1)=sum(sum_Var3(6:13,1))/sum(station_counter(6:13,1)); % Asia
C_Ave_Var3(8,1)=sum(sum_Var3(12:13,1))/sum(station_counter(12:13,1)); % India

end

