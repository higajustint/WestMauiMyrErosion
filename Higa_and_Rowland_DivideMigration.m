%% GENERAL SETUP
clear
clc
close all

% Add path to TopoToolbox 2.0 (Schwanghart and Scherler, 2014)
addpath(genpath('PATH'));
% Add path to DivideTools (Forte and Whipple, 2018)
addpath(genpath('PATH'));

%% Input DEM and shapefile

% 10-m EDNA USGS digital elevation model cropped to West Maui Volcano
dem = GRIDobj('usgs_dem_10m_Wmaui_UTM.asc'); %m

% Set analyses extent
Xex = [739324.203150000,763034.203150000]; % Easting extent [m]
Yex = [2298638.92868800,2327918.92868800]; % Northing extent [m]

% West Maui Volcano caldera boundary adapted from Sherrod et al. (2021)
wai_c_full = shaperead('.\shapefiles\wailuku_caldera_full.shp');

%% Basin outlet point inputs; same as basin maker

names0 = ["Honokohau","Waihee","Iao","Waikapu","Ukumehame",...
    "Olowalu","Launiupoko","Kauaula"]; % In order of below coordinates

% For divide migration at alluvial interface
x0 = [7.48594,7.58002,7.58541,7.582,7.51051,7.485,7.468,7.46]*1e5;
y0 = [2.3262,2.317824,2.311329,2.308244,2.302938,2.3044,2.307535,2.31007]*1e6;

% For divide migration at same elevation (275 m here)
chi_elev = 275; % m
x1 = [7.50854,7.5453,7.55508,7.5725,7.52395,7.50745,7.471415,7.455785]*1e5;
y1 = [2.31944,2.31579,2.310923,2.308275,2.305395,2.30682,2.307636,2.3100035]*1e6;

% For divide migration at same elevation (50 m here)
% chi_elev = 50; % m
% x1 = [7.489578,7.577905,7.60065,7.621727,7.51188,7.48336,7.452022,7.434452]*1e5;
% y1 = [2.324096,2.3178333,2.312611,2.3082235,2.303085,2.304166,2.306296,2.30905]*1e6;

%% Variables

theta = 0.52; % Concavity [unitless]
ref_drain = 10000; % Reference drainage area [m^2]
r_rad = 100; % Relief radius [m]

%%%%%%%%%%%%%%%%%Adjustable variables above only%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create STREAMobj using outlets above

% Crop dem to prescribed values
dem = crop(dem,Xex,Yex); % Crop dem to prescribed extent

FD  = FLOWobj(dem);
ST1 = STREAMobj(FD,'minarea',1000,'unit','mapunits'); % NOT FOR ANALYSES, just to snap points to stream outlets

[~,~,IX0] = snap2stream(ST1,x0,y0,'plot',false); % For bedrock alluvial divide
[~,~,IX1] = snap2stream(ST1,x1,y1,'plot',false); % For same elevation outlets

%% Create drainage basins

DB = drainagebasins(FD,IX0); DB.Z = double(DB.Z); DB.Z(DB.Z == 0) = NaN;

[x00,y00] = ind2coord(dem,IX0); % For bedrock alluvial divide
[x01,y01] = ind2coord(dem,IX1); % For same elevation outlets

%% Order basins for analyses

% compare basin by index in x0 and y0, b1 and b2 are same length
b1 = [1,1,2,3,6,8,5,6,7]; % basin 1 index in x0 and y0
b2 = [2,3,3,4,3,3,4,5,6]; % basin 2 index in x0 and y0

dat1 = cell([length(b1),36]);
dat2 = cell([length(b1),12]);

%% Loop AcrossDivideJH

%%%% AcrossDivide using std_dev wl_method with different elevations %%%%%%%%%%%%

[AREA_OUT]=DivideStability(dem,FD,'verbose',true,'ref_area',ref_drain,'rlf_rad',r_rad,'outlet_control','elevation','min_elevation',0,'theta_ref',theta);

f = waitbar(0,'Please wait...');
for i = 1:length(b1)
    
    % Create input outlets for AcrossDivide function
    rm0 = [x00(b1(i)),y00(b1(i)),1;... % bedrock alluvial
        x00(b2(i)),y00(b2(i)),2];
    
    % AcrossDivide using std_dev wl_method
    [head_vals,out_str_s1]=AcrossDivideJH(dem,FD,AREA_OUT,'outlet_method',...
        'picked_outlets','river_mouths',rm0,'divide_buffer','ultralax',...
        'wl_method','std_dev','display_map',false);close;
    
    % Hard coded to compare only 2 valleys at a time, index valley heads
    % for each side of the divide
    ix1 = head_vals(:,7) == 1;
    ix2 = head_vals(:,7) == 2;
    
    % Mean and standard deviation of elevation on either side of divide
    mean_e1 = mean(head_vals(ix1,3)); s1_e1 = std(head_vals(ix1,3));
    mean_e2 = mean(head_vals(ix2,3)); s1_e2 = std(head_vals(ix2,3));
    
    % Mean and standard deviation of relief on either side of divide
    mean_r1 = mean(head_vals(ix1,5)); s1_r1 = std(head_vals(ix1,5));
    mean_r2 = mean(head_vals(ix2,5)); s1_r2 = std(head_vals(ix2,5));
    
    % Mean and standard deviation of slope on either side of divide
    mean_g1 = mean(head_vals(ix1,4)); s1_g1 = std(head_vals(ix1,4));
    mean_g2 = mean(head_vals(ix2,4)); s1_g2 = std(head_vals(ix2,4));
    
    % Mean and standard deviation of chi on either side of divide
    mean_c1 = mean(head_vals(ix1,6)); s1_c1 = std(head_vals(ix1,6));
    mean_c2 = mean(head_vals(ix2,6)); s1_c2 = std(head_vals(ix2,6));
    
    % Rerun AcrossDivide using ttest wl_method instead of std_dev
    [head_vals_tt,out_str_tt]=AcrossDivideJH(dem,FD,AREA_OUT,'outlet_method',...
        'picked_outlets','river_mouths',rm0,'divide_buffer','ultralax',...
        'wl_method','ttest','display_map',false);close;

    % Hard coded to compare only 2 valleys at a time, index valley heads
    % for each side of the divide (DO AGAIN for new AcrossDivide run)
    ix1 = head_vals_tt(:,7) == 1;
    ix2 = head_vals_tt(:,7) == 2;
    
    % Two sample t-tests against elevation, slope, relief, and chi
    [he,pe] = ttest2(head_vals_tt(ix1,3),head_vals_tt(ix2,3),'Vartype','unequal');
    [hr,pr] = ttest2(head_vals_tt(ix1,5),head_vals_tt(ix2,5),'Vartype','unequal');
    [hg,pg] = ttest2(head_vals_tt(ix1,4),head_vals_tt(ix2,4),'Vartype','unequal');
    [hc,pc] = ttest2(head_vals_tt(ix1,6),head_vals_tt(ix2,6),'Vartype','unequal');
    
    % Make data table1
    dat1_1 = {names0(b1(i)),names0(b2(i)),sum(ix1),sum(ix2),...
        mean_e1,s1_e1,mean_e2,s1_e2,out_str_s1{1},...
        out_str_tt{1},...
        he,pe,...
        mean_r1,s1_r1,mean_r2,s1_r2,out_str_s1{2},...
        out_str_tt{2},...
        hr,pr,...
        mean_g1,s1_g1,mean_g2,s1_g2,out_str_s1{3},...
        out_str_tt{3},...
        hg,pg,...
        mean_c1,s1_c1,mean_c2,s1_c2,out_str_s1{4},...
        out_str_tt{4},...
        hc,pc};
    for j = 1:length(dat1_1)
        dat1{i,j} = dat1_1{j};
    end

    waitbar(i/length(b1),f,'Filling first table');
end
close(f);
%%
%%%% AcrossDivide using std_dev wl_method with same elevation %%%%%%%%%%%%

[AREA_OUT]=DivideStability(dem,FD,'verbose',true,'ref_area',ref_drain,'rlf_rad',r_rad,'outlet_control','elevation','min_elevation',chi_elev,'theta_ref',theta);

f = waitbar(0,'Please wait...');
for i = 1:length(b1)

    % Create input outlets for AcrossDivide function
    rm1 = [x01(b1(i)),y01(b1(i)),1;... % same elevation
        x01(b2(i)),y01(b2(i)),2];

    b1_elev = dem.Z(IX1(b1(i))); % Check that elevations of valley heads
    b2_elev = dem.Z(IX1(b2(i))); % are the same or similar

    [head_vals_1,out_str_s1_1]=AcrossDivideJH(dem,FD,AREA_OUT,'outlet_method',...
        'picked_outlets','river_mouths',rm1,'divide_buffer','ultralax',...
        'wl_method','std_dev','display_map',false); close;

    % Hard coded to compare only 2 valleys at a time, index valley heads
    % for each side of the divide (DO AGAIN for new AcrossDivide run)
    ix1_1 = head_vals_1(:,7) == 1;
    ix2_1 = head_vals_1(:,7) == 2;

    % Mean and standard deviation of chi on either side of divide
    mean_c1_1 = mean(head_vals_1(ix1_1,6)); s1_c1_1 = std(head_vals_1(ix1_1,6));
    mean_c2_1 = mean(head_vals_1(ix2_1,6)); s1_c2_1 = std(head_vals_1(ix2_1,6));

    % Rerun AcrossDivide using ttest wl_method instead of std_dev
    [head_vals_tt_1,out_str_tt_1]=AcrossDivideJH(dem,FD,AREA_OUT,'outlet_method',...
        'picked_outlets','river_mouths',rm1,'divide_buffer','ultralax',...
        'wl_method','ttest','display_map',false); close;

    % Hard coded to compare only 2 valleys at a time, index valley heads
    % for each side of the divide (DO AGAIN for new AcrossDivide run)
    ix1_1 = head_vals_tt_1(:,7) == 1;
    ix2_1 = head_vals_tt_1(:,7) == 2;

    [hc_1,pc_1] = ttest2(head_vals_tt_1(ix1_1,6),head_vals_tt_1(ix2_1,6),'Vartype','unequal');
    
    % Make data table2
    dat2_1 = {b1_elev,b2_elev,sum(ix1_1),sum(ix2_1),...
        mean_c1_1,s1_c1_1,mean_c2_1,s1_c2_1,out_str_s1_1{4},...
        out_str_tt_1{4},...
        hc_1,pc_1};
    for j = 1:length(dat2_1)
        dat2{i,j} = dat2_1{j};
    end

    waitbar(i/length(b1),f,'Filling second table');
end
close(f);
%% Compile data table

dat_tab = [dat1,dat2];
tab = cell2table(dat_tab,"VariableNames",{'Basin 1','Basin 2','n1','n2'...
    'Mean channel head elevation 1 [m]', 'Standard deviation elevation 1 [m]',...
    'Mean channel head elevation 2 [m]', 'Standard deviation elevation 2 [m]','Elevation (standard deviation method)',...
    'Elevation (t-test method)','Elevation t-test h-value','Elevation t-test p-value',...
    'Mean channel head relief 1 [m]', 'Standard deviation relief 1 [m]',...
    'Mean channel head relief 2 [m]', 'Standard deviation relief 2 [m]','Relief (standard deviation method)',...
    'Relief (t-test method)','Relief t-test h-value','Relief t-test p-value',...
    'Mean channel head slope 1 [m/m]', 'Standard deviation slope 1 [m/m]',...
    'Mean channel head slope 2 [m/m]', 'Standard deviation slope 2 [m/m]','slope (standard deviation method)',...
    'slope (t-test method)','slope t-test h-value','slope t-test p-value',...
    'Mean channel head chi 1 [m]', 'Standard deviation chi 1 [m]',...
    'Mean channel head chi 2 [m]', 'Standard deviation chi 2 [m]','Chi (standard deviation method)',...
    'Chi (t-test method)','Chi t-test h-value','Chi t-test p-value',...
    'Basin 1 outlet elevation [m]','Basin 2 outlet elevation [m]',...
    'n1 chi by elevation','n2 chi by elevation'...
    'Mean channel head chi by elevation 1 [m]', 'Standard deviation chi by elevation 1 [m]',...
    'Mean channel head chi by elevation 2 [m]', 'Standard deviation chi by elevation 2 [m]','Chi by elevation (standard deviation method)',...
    'Chi by elevation (t-test method)','Chi by elevation t-test h-value','Chi by elevation t-test p-value',});

%% Make streams for DIVIDEobj
DEM1 = dem; DEM1 = crop(DEM1,[739324.29057,758500],[2298639.86446800,2327918.84126800]); % Crop to exclude alluvial areas for DIVIDEobj
FD1  = FLOWobj(DEM1);

A = flowacc(FD1).*(FD1.cellsize^2);

ST2 = STREAMobj(FD1,A>ref_drain); % create streamobj for channels > ref_drain

%% Create clean DIVIDEobj

D = DIVIDEobj(FD1,ST2);
D2 = cleanedges(D,FD1);
D = divorder(D2,'strahler');

%% Shrink divide object by order

Dtest = shrink(D,FD1,'order',4);

%% Across divide difference in hillslope relief (m)

DZ = vertdistance2stream(FD1,ST2,DEM1);
DZ.Z(isinf(DZ.Z)) = nan;

%% Calculate divide assymetry index

[MS,S] = asymmetry(Dtest,DZ);
for i = 1 : length(S)
  S(i).length = max(getdistance(S(i).x,S(i).y));
end

%% Plot divide assymetry index with arrows

fig = figure(1);
tiledlayout(2,4,"TileSpacing","compact");
nexttile([2,2])
imageschs(dem,[],'colormap',[.9 .9 .9],'colorbar',false,'ticklabels','nice');
hold on
plotc(Dtest,vertcat(S.rho),'caxis',[0 0.5],'limit',[1000 inf])
colormap(gca,flipud(pink(5)))
axis image
hc = colorbar;
hc.Label.String = 'Divide asymmetry index';
mapshow(wai_c_full,'LineStyle',':','FaceColor','none','EdgeColor','k','LineWidth',2);
plot([749020.911288448,751209.075065469,751209.075065469,749020.911288448,749020.911288448],[2307020.27548124,2307020.27548124,2309722.60787060,2309722.60787060,2307020.27548124],'k-','LineWidth',1);
plot([749722.719966948,753599.185575932,753599.185575932,749722.719966948,749722.719966948],[2310807.31903960,2310807.31903960,2315594.66570963,2315594.66570963,2310807.31903960],'k-','LineWidth',1);
plot([746941.243505339,757260.394956453,757260.394956453,746941.243505339,746941.243505339],[2304836.83634667,2304836.83634667,2312208.79476447,2312208.79476447,2304836.83634667],'k-','LineWidth',1);
set(gca,'XTick',749484);
set(gca,'XTickLabel',['156.6' char(176) 'W']);
set(gca,'YTick',[2302030,2324182]);
set(gca,'YTickLabel',[['20.8' char(176) 'N'];['21.0' char(176) 'N']]);
plot([759384,762384],[2299282,2299282],'k-','LineWidth',3);
text(760134,2300000,'3 km');
patch([759384,760884,762384,760884],[2323182,2327182,2323182,2324182],rgb('black'),'EdgeColor','none');
title('Drainage divide asymmetry')

nexttile;
imageschs(dem,[],'colormap',[.9 .9 .9],'colorbar',false,'ticklabels','none');
hold on
plotc(Dtest,vertcat(S.rho),'caxis',[0 0.5],'limit',[1000 inf])
plot([750091,751091],[2307180,2307180],'k-','LineWidth',2);
colormap(gca,flipud(pink(5)))
axis image
ix = [MS.order]>1 & [S.length]>1000; % Plot red arrows for divides > order 1 and length 1000 m
quiver([MS(ix).X],[MS(ix).Y],[MS(ix).u],[MS(ix).v],1,...
  'color','r','linewidth',2)
xlim([749020.911288448,751209.075065469]);
ylim([2307020.27548124,2309722.60787060]);

nexttile;
imageschs(dem,[],'colormap',[.9 .9 .9],'colorbar',false,'ticklabels','none');
hold on
plotc(Dtest,vertcat(S.rho),'caxis',[0 0.5],'limit',[1000 inf])
plot([752353,753353],[2311044,2311044],'k-','LineWidth',2);
colormap(gca,flipud(pink(5)))
axis image
ix = [MS.order]>1 & [S.length]>1000; % Plot red arrows for divides > order 1 and length 1000 m
quiver([MS(ix).X],[MS(ix).Y],[MS(ix).u],[MS(ix).v],1,...
  'color','r','linewidth',2)
xlim([749722.719966948,753599.185575932]);
ylim([2310807.31903960,2315594.66570963]);

nexttile([1,2]);
imageschs(dem,[],'colormap',[.9 .9 .9],'colorbar',false,'ticklabels','none');
hold on
plotc(Dtest,vertcat(S.rho),'caxis',[0 0.5],'limit',[1000 inf])
plot([756192,757192],[2305371,2305371],'k-','LineWidth',2);
colormap(gca,flipud(pink(5)))
axis image
ix = [MS.order]>1 & [S.length]>1000; % Plot red arrows for divides > order 1 and length 1000 m
mapshow(wai_c_full,'LineStyle',':','FaceColor','none','EdgeColor','k','LineWidth',2);
quiver([MS(ix).X],[MS(ix).Y],[MS(ix).u],[MS(ix).v],1.5,...
  'color','r','linewidth',2)
xlim([746941.243505339,757260.394956453]);
ylim([2304836.83634667,2312208.79476447]);

fig.Position = [212.5555555555555,306.3333333333333,1157.333333333333,545.7777777777778];