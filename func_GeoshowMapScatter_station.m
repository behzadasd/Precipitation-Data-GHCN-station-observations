function [h1]=func_GeoshowMapScatter_station(Variable, Lat, Lon, var_limit, x_limit, y_limit, dir_out_name_fig, title_text)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load coast
geoshow(flipud(lat), flipud(long),'DisplayType', 'polygon', 'FaceColor', [0.94 0.97 1])
%geoshow('landareas.shp', 'FaceColor', [0.97 0.97 0.97]);  %%% Enable this to Change the LandArea color
states = shaperead('usastatehi', 'UseGeoCoords', true);    %%% Enable this to add USA state boundaries
geoshow(states,'DefaultFaceColor', 'white', 'DefaultEdgeColor', 'blue', 'FaceColor', 'white');
xlim(x_limit); ylim(y_limit);

hold on

h1=scatter(Lon, Lat, ones(size(Lon,1),1)*10 , Variable , 'fill');
%set(h1,'alphadata',~isnan(Variable)) % Sets NaN values no color (colors them white)
set(gca, 'CLim', var_limit);
colormap jet
title(title_text)
set(gca,'YDir','normal') % Prevents fliping latitude directions because of the upper arrays being greater than lower arrays
colorbar('location','Eastoutside')
set(gca,'FontSize',34, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
saveas(h1,dir_out_name_fig)
saveas(h1,dir_out_name_fig,'png')
%close


end
