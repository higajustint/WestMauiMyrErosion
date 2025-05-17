%% GENERAL SETUP
clear
clc
close all

% Add path to TopoToolbox 2.0 (Schwanghart and Scherler, 2014)
addpath(genpath('PATH'));
% Add path to RGB triple of color name, version 2 (Jonasson, 2025)
addpath(genpath('PATH'));

%% Inputs

% 10-m EDNA USGS digital elevation model cropped to West Maui Volcano
dem = GRIDobj('FILE'); %m

% Distance raster from Bouguer gravity high (Strange et al., 1965)
dist2cal = GRIDobj('Strange_grav_cald_dist.asc'); % m

% Hawaii Rainfall Atlas and standard deviation (Giambelluca et al., 2013)
rainatlas = GRIDobj('FILE'); % mm/yr
rainatlas_std = GRIDobj('FILE'); % mm/yr

% Valley ID number, width [m], and length [m] table
v_wid = readtable('Maui_valley_aspect.csv');

% Relict surface from thin-plate spline
surface1 = GRIDobj('Wmaui_thinplatespline.asc'); % m

% Set analyses extent
Xex = [739324.203150000,763034.203150000]; % Easting extent [m]
Yex = [2298638.92868800,2327918.92868800]; % Northing extent [m]

%% Shapefiles

% All basins outlined
basins = shaperead('.\shapefiles\WMV_Valleys.shp');

% Basins named in paper
basins_MA = shaperead('.\shapefiles\WMV_Valleys_main.shp');

% West Maui Volcano caldera boundary adapted from Sherrod et al. (2021)
wai_c_full = shaperead('.\shapefiles\Wailuku_caldera_full.shp');

%% Variables

% Want to regress only valleys larger than a certain size?
size_thresh = 0; % km2 % NOTE if not 0, fprintf section will fail, proceed anyway

% Incision threshold for streams
area_thresh = 1e5; % m^2 for streams

% Resample resolution
resampres = 10; % m

% Radius for relief calculation
radius = 1000; % m

% Average density of rock
rho = 3000; % m^3/kg from Ferrier et al. (2013, GSA Bulletin)
rho_km = rho*1e9; % kg/m^3 to kg/km^3

% What type of rainfall statistic do you want? "Mean", "Median", "Min", "Max", "percentile"?
% If percentile, fill in rainprc = percentile value
rainstr = 'Mean'; rainprc = NaN; % set as mean annual precipitation MAP

% What type of distance to caldera statistic do you want? "Mean", "Median", "Min", "Max", "percentile"?
% If percentile, fill in distprc = percentile value
diststr = 'Min'; distprc = NaN; % set as minimum distance to volcanic center Dmin

% What type of relief statistic do you want? "Mean", "Median", "Min", "Max", "percentile"?
% If percentile, fill in relprc = percentile value
relstr = 'Mean'; relprc = NaN; % set as mean relief

% What type of slope statistic do you want? "Mean", "Median", "Min", "Max", "percentile"?
% If percentile, fill in slpprc = percentile value
slpstr = 'Mean'; slpprc = NaN; % set as mean slope

% What type of Relict Slope statistic do you want? "Mean", "Median", "Min", "Max", "percentile"?
% If percentile, fill in slpprc = percentile value
s0str = 'percentile'; s0prc = 95; % set as 95th steepest percentile of relict slope

% What type of valley width statistic do you want? "Mean", "Median", "Min", "Max", "percentile"?
% If percentile, fill in slpprc = percentile value
widstr = 'Max'; widprc = NaN; % set as max valley width

% Standardize PCA? Yes or No?
pca_stand = 'Yes';

% Color PCA by "distance 2 caldera", "precipitation", "erosion percentile"?
% If "erosion percentile", do you want er_type "valley volume", "fraction eroded", or 
% "erosion rate"?
pca_col = 'distance 2 caldera'; er_type = NaN;

% PCA variables

% For precipitation or distance to caldera
varlab = {'R','S','Asp','A'}; % names [relief,slope,aspect ratio, area]
arr_vars = [7,8,11,6]; % arr matrix column index numbers that match names %% NOT SAME AS TAB COLUMN INDEX

%% Basin outlet point inputs; same as basin maker

% Basin names
Bname = ["Kauaula tributary 1","Kanaha","Kahoma","Wahikuli","Honokowai","Kahana","Honokahua","Honolua","Honokohau no relict","Anakaluahine","Poelua","Hononana","Waihali","Kahakuloa","Waipili","Wailena","Makamakaole","Waihee","Kalaeiliili","Kalepa","Kope 1","Kope 2","North Waiehu","South Waiehu","Waiehu","Iao","Waikapu","Waikapu tributary 1","Waikapu tributary 2","Pohakea","Papalaua","Ukumehame","Olowalu","Launiupoko","Launiupoko tributary 1","Launiupoko tributary 2","Luakoi","Kauaula"];

% outlet of basins
x0 = [745892.410030000	742954.687190000	743003.812990000	740449.271390000	741490.738350000	742198.149870000	744182.832190000	745489.578470000	748594.329030000	750775.514550000	751591.002830000	753025.476190000	754440.299230000	754803.830150000	756071.275790000	756208.828030000	757299.420790000	757997.007150000	758124.734230000	758399.838710000	758193.510350000	757997.007150000	758095.258750000	757928.231030000	757987.181990000	758537.390950000	758222.985830000	757800.503950000	758105.083910000	758733.894150000	752966.525230000	751050.619030000	748496.077430000	746815.975070000	746452.444150000	746098.738390000	745460.102990000	746000.486790000];
y0 = [2310248.29100800	2311869.44240800	2312154.37204800	2314306.08208800	2317961.04160800	2321183.69408800	2324072.29112800	2325398.68772800	2326204.35084800	2326489.28048800	2325909.59604800	2325369.21224800	2323993.68984800	2323797.18664800	2322559.21648800	2322441.31456800	2320211.00324800	2317823.48936800	2317135.72816800	2316673.94564800	2315249.29744800	2314915.24200800	2314080.10340800	2313726.39764800	2313215.48932800	2311329.05860800	2308253.78352800	2306495.07988800	2305846.61932800	2303577.00736800	2301700.40180800	2302948.19712800	2304402.32080800	2307546.37200800	2307693.74940800	2308204.65772800	2309118.39760800	2310071.43812800];

% To remove Honokohau relict surface, outlets of relict topography for Honokohau
xo0 = [751359.203150000 750909.203150000];
yo0 = [2313593.92868800 2313623.92868800];

% Indicies for specific valleys
Honokohau = 9;
Waihee = 18;
Iao = 26;
Ukumehame = 32;
Olowalu = 33;

%% Input age of topmost lava flow

% Rough estimate of end of West Maui Volcano major activity 1.1 - 1.3 Ma
age = 1.20; % Ma
age_err = 0.1; % Ma 1sigma

%%%%%%%%%%%%%%%%%Adjustable variables above only%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resample DEM to desired resolution

dem = resample(dem,resampres,'bicubic');
dem = crop(dem,Xex,Yex); % Crop dem to prescribed extent
[dxd,dyd] = refmat2XY(dem.refmat, dem.size); % Convert referencing matrix (refmat) to coordinate vectors
[dxg,dyg] = meshgrid(dxd,dyd);

%% Interpolate distance to caldera

[d2xd,d2yd] = refmat2XY(dist2cal.refmat, dist2cal.size); % Convert referencing matrix (refmat) to coordinate vectors
[d2xg,d2yg] = meshgrid(d2xd,d2yd);

dist2calZ = interp2(d2xg,d2yg,dist2cal.Z,dxg,dyg,'cubic');
dist2cal = dem; dist2cal.Z = dist2calZ;

% coordinate of closest point
[xxmin,yymin] = ind2coord(dist2cal,find(dist2cal.Z == min(dist2cal.Z(:))));

%% Interpolate rainfall atlas and standard deviation

[rxd,ryd] = refmat2XY(rainatlas.refmat, rainatlas.size); % Convert referencing matrix (refmat) to coordinate vectors
[rxg,ryg] = meshgrid(rxd,ryd);

rainZ = interp2(rxg,ryg,rainatlas.Z,dxg,dyg,'cubic');
rainatlas = dem; rainatlas.Z = rainZ;

[rsxd,rsyd] = refmat2XY(rainatlas_std.refmat, rainatlas_std.size); % Convert referencing matrix (refmat) to coordinate vectors
[rsxg,rsyg] = meshgrid(rsxd,rsyd);

rainstdZ = interp2(rsxg,rsyg,rainatlas_std.Z,dxg,dyg,'cubic');
rainatlas_std = dem; rainatlas_std.Z = rainstdZ;

%% Interpolate minimum and maximum reconstructed surfaces

[m1xd,m1yd] = refmat2XY(surface1.refmat, surface1.size); % Convert referencing matrix (refmat) to coordinate vectors
[m1xg,m1yg] = meshgrid(m1xd,m1yd);

surfZ = interp2(m1xg,m1yg,surface1.Z,dxg,dyg,'cubic');
surface1 = dem; surface1.Z = surfZ;

%% Difference dem minus relict

diff_surf1 = dem; diff_surf1.Z = surface1.Z - dem.Z;

%% Calculate basin relief and slope

relf = localtopography(dem,radius,'type','range'); % m

slp = arcslope(dem,'deg');

%% Calculate slope of original surface (i.e., max surface)
s0 = arcslope(surface1,'deg');
s0 = localtopography(s0,radius,'type','mean');
E0 = localtopography(surface1,radius,'type','mean');

%% Calculate streams

FD = FLOWobj(dem);
DIST = flowdistance(FD);
A = flowacc(FD).*(FD.cellsize^2);
S = STREAMobj(FD,A>area_thresh);

[xsnap,ysnap,IX] = snap2stream(S,x0,y0,'plot',true);

DB = drainagebasins(FD,IX); DB.Z = double(DB.Z); DB.Z(DB.Z == 0) = NaN;

if ~isempty(xo0)
    [xsnapo0,ysnapo0,IXo0] = snap2stream(S,xo0,yo0,'plot',true);
    DBo0 = drainagebasins(FD,IXo0); DBo0.Z = double(DBo0.Z); DBo0.Z(DBo0.Z == 0) = NaN;
    % Make all Honokohau relict topo NaN
    DB.Z(~isnan(DBo0.Z)) = NaN;
end

%% Basin averaged Ksn from TopoToolbox website
S1 = STREAMobj(FD,'outlets',IX,'minarea',area_thresh,'unit','mapunits');

% Concavity theta (m/n)
mn_optimized = mnoptim(S1,dem,flowacc(FD),'mnrange',[0.1 1],'crossval',false,...
    'plot',false,'a0',area_thresh);mn_optimized = mn_optimized{1,1};
close all;

m = mn_optimized; % Assuming n = 1 (e.g., Ferrier et al., 2013, Nature)
m = round(m,2); % Round m to the nearest hundredth
n = 1;
theta = m/n; % Use optimized theta

d1   = getnal(S1,DB); d = d1; d(isnan(d)) = [];
k   = ksn(S1,dem,flowacc(FD),theta); k(isnan(d1)) = [];
km  = accumarray(d,k,[],@mean);
km_std  = accumarray(d,k,[],@std);
K   = GRIDobj(dem)*nan;
K.Z(DB.Z>0) = km(DB.Z(DB.Z>0));
K_std   = GRIDobj(dem)*nan;
K_std.Z(DB.Z>0) = km_std(DB.Z(DB.Z>0));

%% KsnQ

zb1=fillmissing(rainatlas.Z, 'linear');         % Fillmissing() substitutes NaNs partially
zb1(zb1<0) = 0;
raintest = rainatlas; raintest.Z = zb1./1000; % m

% Turn rainfall into discharge
dischm = flowacc(FD,raintest.Z).*(FD.cellsize).^2; %m3/yr

%
D8 = gradient8(dem,'tan');
KSNQ = D8.Z./((dischm).^-theta); % Adams et al. (2020) Ksn-q

% Get the same index for stream pixels as getnal function above
DBidx = ~isnan(DB.Z);
Aidx = A.Z>area_thresh;
DBAidx  = logical(DBidx.*Aidx);

% Make new d vector of DB values
DB_d = DB.Z(DBAidx);

% Filter out all non-stream pixels
KSNQ_thresh = KSNQ.Z(DBAidx);

% Accumarray mean Ksn-q and standard deviation
kQm  = accumarray(DB_d,KSNQ_thresh,[],@(x)mean(x,'omitnan'));
kQm_std  = accumarray(DB_d,KSNQ_thresh,[],@(x)std(x,'omitnan'));

% Apply mean ksnq and standard dev to a GRIDobj
KQ   = GRIDobj(dem)*nan;
KQ.Z(DB.Z>0) = kQm(DB.Z(DB.Z>0));
KQ_std   = GRIDobj(dem)*nan;
KQ_std.Z(DB.Z>0) = kQm_std(DB.Z(DB.Z>0));

%% Make table of all characteristics for basins

arr = NaN(max(DB.Z(:)),35); % Change second number as more columns are needed
if ~isempty(xo0) % If you are calculating without Honokohau relict you need one less row
    relf_val = cell([1,max(DB.Z(:))-1]);
    slp_val = cell([1,max(DB.Z(:))-1]);
    asp_val = cell([1,max(DB.Z(:))-1]);
else
    relf_val = cell([1,max(DB.Z(:))]);
    slp_val = cell([1,max(DB.Z(:))]);
    asp_val = cell([1,max(DB.Z(:))]);
end

counter9 = 0; % Count number of times i forcibly set to Honokohau value to take care of duplicate for relict analyses

f = waitbar(0,'Please wait...');

for i = 1:max(DB.Z(:))

    % Calculate precipitation factor
    rain_idx = rainatlas.Z(DB.Z == i);
    
    if strcmp(rainstr,'Min')
        rain_fac = min(rain_idx(:)); % m
    elseif strcmp(rainstr,'Max')
        rain_fac = max(rain_idx(:)); % m
    elseif strcmp(rainstr,'Mean')
        % Mean
        rain_fac = mean(rain_idx,'omitnan'); % mm/yr
    elseif strcmp(rainstr,'Median')
        % Median
        rain_fac = median(rain_idx,'omitnan'); % mm/yr
    elseif strcmp(rainstr,'percentile')
        % Rainfall xxth percentile
        rain_fac = prctile(rain_idx,rainprc); % mm/yr  
    else
        disp('Must specify Mean, Median, percentile, case sensitive');
        return
    end

    % Calculate eroded volume
    vol_eroded_surf1 = sum(diff_surf1.Z(DB.Z == i),'omitnan')...
        *diff_surf1.cellsize^2/1e9; % m2 to km3

    % Calculate relict "wedge" volume
    relict_vol_surface1 = sum(surface1.Z(DB.Z == i),'omitnan')...
        *surface1.cellsize^2/1e9; % m2 to km3

    % Calculate fraction of wedge eroded
    fx = vol_eroded_surf1/relict_vol_surface1; % km3/km3

    % Calculate basin area
    num_cells = DB.Z == i; num_cells = sum(num_cells(:)); % no. cells
    Ab = num_cells*DB.cellsize^2; % m2
    Abkm2 = Ab/1e6; % m2 to km2

    % Calculate million-year erosion rate
    Ev = rho_km*vol_eroded_surf1/Abkm2/age/(1e3*1e6); % kg/km2/Ma to t/km2/yr
    % Propagage error to get 1sigma
    Ev_1s = age_err/age*Ev; % Ma/Ma * t/km2/yr

    % Calculate vertical erosion rate and error in mm/yr
    Ev_v = Ev/rho; % t/km2/yr*kg/m3 = mm/yr

    Ev_v_1s = Ev_1s/rho;

    % Calculate fraction of vol gained in error after erosion
    diff_vec = diff_surf1.Z(DB.Z == i);
    pos_diff = find(diff_vec > 0); % no. volume gained cells
    neg_diff = find(diff_vec < 0); % no. volume eroded cells

    pos_vol = sum(diff_vec(pos_diff))*diff_surf1.cellsize^2; % m2
    neg_vol = sum(diff_vec(neg_diff))*diff_surf1.cellsize^2; % m2

    vol_err_frac = abs(neg_vol)/pos_vol; % m2/m2

    % Calculate mean, median, or min distance to caldera

    if strcmp(diststr,'Min')
        cal_dist = min(dist2cal.Z(DB.Z == i)); % m
    elseif strcmp(diststr,'Max')
        cal_dist = max(dist2cal.Z(DB.Z == i)); % m
    elseif strcmp(diststr,'Mean')
        cal_dist = mean(dist2cal.Z(DB.Z == i),'omitnan'); % m
    elseif strcmp(diststr,'Median')
        cal_dist = median(dist2cal.Z(DB.Z == i),'omitnan'); % m
    elseif strcmp(diststr,'percentile')
        % Distance > xth percentile
        cal_dist = prctile(dist2cal.Z(DB.Z == i),distprc); % m
    else
        disp('Must specify percentile, Mean, Median, Min, or Max, case sensitive');
        return
    end

    % Calculate mean, median, or min basin relief

    if strcmp(relstr,'Min')
        rel_fac = min(relf.Z(DB.Z == i)); % m
    elseif strcmp(relstr,'Max')
        rel_fac = max(relf.Z(DB.Z == i)); % m
    elseif strcmp(relstr,'Mean')
        rel_fac = mean(relf.Z(DB.Z == i),'omitnan'); % m
    elseif strcmp(relstr,'Median')
        rel_fac = median(relf.Z(DB.Z == i),'omitnan'); % m
    elseif strcmp(relstr,'percentile')
        % Mean of distance > xth percentile
        rel_fac = prctile(relf.Z(DB.Z == i),relprc); % m
    else
        disp('Must specify percentile, Mean, Median, Min, or Max, case sensitive');
        return
    end
    
    % Hold all relief cells within a given valley
    relf_val{i} = relf.Z(DB.Z == i);

    % Calculate mean, median, or min basin slope

    if strcmp(slpstr,'Min')
        slp_fac = min(slp.Z(DB.Z == i)); % deg.
    elseif strcmp(slpstr,'Max')
        slp_fac = max(slp.Z(DB.Z == i)); % deg.
    elseif strcmp(slpstr,'Mean')
        slp_fac = mean(slp.Z(DB.Z == i),'omitnan'); % deg.
    elseif strcmp(slpstr,'Median')
        slp_fac = median(slp.Z(DB.Z == i),'omitnan'); % deg.
    elseif strcmp(slpstr,'percentile')
        % Mean of distance > xth percentile
        slp_fac = prctile(slp.Z(DB.Z == i),slpprc); % deg.
    else
        disp('Must specify percentile, Mean, Median, Min, or Max, case sensitive');
        return
    end

    % Hold all slope cells within a given valley
    slp_val{i} = slp.Z(DB.Z == i);

    % Calculate mean, median, or min original basin slope

    if strcmp(s0str,'Min')
        s0_fac = min(s0.Z(DB.Z == i));
    elseif strcmp(s0str,'Max')
        s0_fac = max(s0.Z(DB.Z == i));
    elseif strcmp(s0str,'Mean')
        s0_fac = mean(s0.Z(DB.Z == i),'omitnan');
    elseif strcmp(s0str,'Median')
        s0_fac = median(s0.Z(DB.Z == i),'omitnan');
    elseif strcmp(s0str,'percentile')
        % Mean of distance > xth percentile
        s0_fac = prctile(s0.Z(DB.Z == i),s0prc);
    else
        disp('Must specify percentile, Mean, Median, Min, or Max, case sensitive');
        return
    end

    % Calculate mean, median, or min basin width
    
    v_idx = find(v_wid{:,1} == i);
    if strcmp(widstr,'Min')
        wid_fac = min(v_wid{:,2}(v_idx)); % m
    elseif strcmp(widstr,'Max')
        wid_fac = max(v_wid{:,2}(v_idx)); % m
    elseif strcmp(widstr,'Mean')
        wid_fac = mean(v_wid{:,2}(v_idx),'omitnan'); % m
    elseif strcmp(widstr,'Median')
        wid_fac = median(v_wid{:,2}(v_idx),'omitnan'); % m
    elseif strcmp(widstr,'percentile')
        % Mean of distance > xth percentile
        wid_fac = prctile(v_wid{:,2}(v_idx),widprc); % m
    else
        disp('Must specify percentile, Mean, Median, Min, or Max, case sensitive');
        return
    end

    % Length of simplified stream as proxy for valley length
    len_fac = mean(v_wid{:,3}(v_idx)); % m all lengths the same, mean is equal to length

    % Valley aspect as a ratio of wid_fac/len_fac
    asp_ratio = wid_fac/len_fac; % m/m

    % Hold all aspect values within a given valley
    asp_val{i} = v_wid{:,2}(v_idx)/len_fac;

    % Index basin averaged ksn and standard deviation in m^(2theta)
    BA_ksn = unique(K.Z(DB.Z == i));
    BA_ksn_1s = unique(K_std.Z(DB.Z == i));

    % Index basin averaged ksn and standard deviation in m^(2theta)
    BA_ksnq = unique(KQ.Z(DB.Z == i));
    BA_ksnq_1s = unique(KQ_std.Z(DB.Z == i));

    % ksn from chi basin averaged m^(2theta)
    S2 = STREAMobj(FD,'outlets',IX(i),'minarea',area_thresh,'unit','mapunits');
    C = chiplot(S2,dem,flowacc(FD),'mn',theta,'a0',1,'plot',false);
    Cbeta_1s = sqrt(sum(d == i))*C.betase; % convert to standard deviation adapted from "Mean basin ksn and smoothing madness" TopoToolbox blog (i.e., use n = # of stream segments for ksn by segment)

    % Erodibility in meters^(1-2m)/yr from chi
    K_chi = Ev_v/C.beta^n/1e3; % convert EvL_v from mm/yr to m/yr
    K_chi_1s = K_chi*sqrt((Ev_v_1s/Ev_v)^2 + (Cbeta_1s/C.beta)^2); %%%%%%% AT THE MOMENT, ALL ERROR CALCULATIONS ASSUME n = 1

    % Erodibility meters^(1-2m)/yr from stream segments
    K_BA = Ev_v/BA_ksn^n/1e3; % convert EvL_v from mm/yr to m/yr
    K_BA_1s = K_BA*sqrt((Ev_v_1s/Ev_v)^2 + (BA_ksn_1s/BA_ksn)^2);

    % Erodibility meters^(1-2m)/yr from stream segments Ksn-q
    Kq_BA = Ev_v/BA_ksnq^n/1e3; % convert EvL_v from mm/yr to m/yr
    Kq_BA_1s = Kq_BA*sqrt((Ev_v_1s/Ev_v)^2 + (BA_ksnq_1s/BA_ksnq)^2);
    
    % Klp/Kq from chi method on flow accumulation of rainfall
    CQ = chiplot(S2,dem,flowacc(FD,raintest.Z),'mn',theta,'a0',1,'plot',false);
    CQbeta_1s = sqrt(sum(d == i))*CQ.betase; % convert to standard deviation adapted from "Mean basin ksn and smoothing madness" TopoToolbox blog (i.e., use n = # of stream segments for ksn by segment)

    Klp_chi = Ev_v/CQ.beta^n/1e3; % convert EvL_v from mm/yr to m/yr
    Klp_chi_1s = Klp_chi*sqrt((Ev_v_1s/Ev_v)^2 + (CQbeta_1s/CQ.beta)^2); %%%%%%% AT THE MOMENT, ALL ERROR CALCULATIONS ASSUME n = 1

    % Klp erodibility and standard deviation (Adams et al., 2020)
    Rmean = mean(rain_idx,'omitnan')/1000; % mm/yr to m/yr calculate mean precipitation of Area (may be duplicate)
    Klp_MAP = K_chi/Rmean^m; % m^(1-2m)/yr*(yr/m)^m = m^(1-3m)yr^(m-1)
    Klp_MAP_1s = sqrt( ((1/(Rmean^m))*K_chi)^2 + ((-(K_chi*m/(Rmean^(m+1))))*std(rain_idx,'omitnan'))^2 );

    % Fill array to store data

    arr(i,1) = i;
    arr(i,2) = xsnap(i);
    arr(i,3) = ysnap(i);
    arr(i,4) = cal_dist;
    arr(i,5) = rain_fac;
    arr(i,6) = Abkm2;
    arr(i,7) = rel_fac;
    arr(i,8) = slp_fac;
    arr(i,9) = wid_fac;
    arr(i,10) = len_fac;
    arr(i,11) = asp_ratio;
    arr(i,12) = vol_eroded_surf1;
    arr(i,13) = vol_err_frac;
    arr(i,14) = fx;
    arr(i,15) = Ev;
    arr(i,16) = Ev_1s;
    arr(i,17) = Ev_v;
    arr(i,18) = Ev_v_1s;
    arr(i,19) = BA_ksn;
    arr(i,20) = BA_ksn_1s;
    arr(i,21) = C.beta;
    arr(i,22) = Cbeta_1s;
    arr(i,23) = K_BA*1e-0;
    arr(i,24) = K_BA_1s*1e-0;
    arr(i,25) = K_chi*1e-0;
    arr(i,26) = K_chi_1s*1e-0;
    arr(i,27) = BA_ksnq;
    arr(i,28) = Kq_BA*1e-0;
    arr(i,29) = Kq_BA_1s*1e-0;
    arr(i,30) = CQ.beta;
    arr(i,31) = Klp_chi*1e-0;
    arr(i,32) = Klp_chi_1s*1e-0;
    arr(i,33) = Klp_MAP*1e-0;
    arr(i,34) = Klp_MAP_1s*1e-0;
    arr(i,35) = s0_fac;

    waitbar(i/max(DB.Z(:)),f,'Filling table');

end
close(f);

% Make table
tab = array2table(arr,"VariableNames",["num","Easting","Northing",...
    "cal_dist","rain_fac",...
    "Abkm2","rel_fac","slp_fac", "wid_fac","len_fac","asp_ratio",...
    "vol_eroded_L","vol_err_frac_L",...
    "fx_L",...
    "EvL","EvL_1s",...
    "EvL_v","EvL_v_1s",...
    "BA_ksn", "BA_ksn_1s",...
    "Cbeta","Cbeta_1s",...
    "KL_BA10e6","KL_BA_1s10e6",...
    "KL_chi10e6","KL_chi_1s10e6",...
    "BAksnq",...
    "KlpLBA10e6","KlpLBA1s10e6",...
    "chiksnq",...
    "KlpLchi10e6","KlpLchi1s10e6",...
    "KlpLMAP10e6","KlpLMAP1s10e6",...
    "OriginalSlope"]);

tab = addvars(tab,Bname','Before','num','NewVariableNames',{'Name'});

%% linear comparison for ksn and K from different methods

% Chi is first variable, stream segments or Klp is second
ksn_fit = fitlm(tab.Cbeta,tab.BA_ksn);

%% Organize valleys by quartiles for PCA plotting

if strcmp(pca_col,'erosion percentile')
    prctileG_max = DB; % Color grid object basins by percentile
    prctileG_max.Z(~isnan(prctileG_max.Z)) = 2; % set interquartile range value = 2

    prctileG_min = DB; % Color grid object basins by percentile
    prctileG_min.Z(~isnan(prctileG_min.Z)) = 2; % set interquartile range value = 2

    if strcmp(er_type,'valley volume')
        top_qrtle_min = prctile(tab.vol_eroded_L,75); % Calculate top quartile

        bot_qrtle_min = prctile(tab.vol_eroded_L,25);

        top_idx_min = arr(tab.vol_eroded_L > top_qrtle_min,1); % Index those > top quartile

        bot_idx_min = arr(tab.vol_eroded_L < bot_qrtle_min,1); % Index those < bottom quartile

    elseif strcmp(er_type,'fraction eroded')
        top_qrtle_min = prctile(tab.fx_L,75);

        bot_qrtle_min = prctile(tab.fx_L,25);

        top_idx_min = arr(tab.fx_L > top_qrtle_min,1);

        bot_idx_min = arr(tab.fx_L < bot_qrtle_min,1);

    elseif strcmp(er_type,'erosion rate')
        top_qrtle_min = prctile(tab.EvL,75);

        bot_qrtle_min = prctile(tab.EvL,25);

        top_idx_min = arr(tab.EvL > top_qrtle_min,1);

        bot_idx_min = arr(tab.EvL < bot_qrtle_min,1);
    end
    
    for i = 1:length(top_idx_min)
        prctileG_min.Z(DB.Z == top_idx_min(i)) = 1;
    end

    for i = 1:length(bot_idx_min)
        prctileG_min.Z(DB.Z == bot_idx_min(i)) = 3;
    end

    for i = 1:length(top_idx_max)
        prctileG_max.Z(DB.Z == top_idx_max(i)) = 1;
    end

    for i = 1:length(bot_idx_max)
        prctileG_max.Z(DB.Z == bot_idx_max(i)) = 3;
    end

elseif strcmp(pca_col,'distance 2 caldera')
    prctileG = DB; % Color grid object basins by percentile
    prctileG.Z(~isnan(prctileG.Z)) = 2; % Set interquartile range value = 2

    top_qrtle = prctile(tab.cal_dist,75); % Calculate top quartile
    bot_qrtle = prctile(tab.cal_dist,25); % Calculate bottom quartile

    top_idx = arr(tab.cal_dist > top_qrtle,1);
    bot_idx = arr(tab.cal_dist < bot_qrtle,1);

    for i = 1:length(top_idx)
        prctileG.Z(DB.Z == top_idx(i)) = 3;
    end

    for i = 1:length(bot_idx)
        prctileG.Z(DB.Z == bot_idx(i)) = 1;
    end

elseif strcmp(pca_col,'precipitation')
    prctileG = DB; % Color grid object basins by percentile
    prctileG.Z(~isnan(prctileG.Z)) = 2; % Set interquartile range value = 2

    top_qrtle = prctile(tab.rain_fac,75); % Calculate top quartile
    bot_qrtle = prctile(tab.rain_fac,25); % Calculate bottom quartile

    top_idx = arr(tab.rain_fac > top_qrtle,1);
    bot_idx = arr(tab.rain_fac < bot_qrtle,1);

    for i = 1:length(top_idx)
        prctileG.Z(DB.Z == top_idx(i)) = 1;
    end

    for i = 1:length(bot_idx)
        prctileG.Z(DB.Z == bot_idx(i)) = 3;
    end
else
    disp('Input PCA color is not an option, case sensitive');
    return
end

%% Make GRIDobj of volume eroded, fraction eroded, and erosion rate

% Dist 2 caldera GRIDobj
D2Cgrid = DB;

% Morphologic metrics GRIDobj
AreaA = DB;
AspectAsp = DB;
ReliefR = DB;
SlopeS = DB;

% Volume eroded GRIDobj
vol_gridmin = DB;

% Fraction eroded GRIDobj
frac_gridmin = DB;

% Erosion rate GRIDobj
rate_gridmin = DB;

% Neg vol error fraction GRIDobj
neg_gridmin = DB;

% K from chi
Kchi_gridmin = DB;

% Klp
Klp_gridmin = DB;

for i = 1:max(DB.Z(:))

    D2Cgrid.Z(DB.Z == i) = tab.cal_dist(i);

    % Surface
    
    vol_gridmin.Z(DB.Z == i) = tab.vol_eroded_L(i);

    frac_gridmin.Z(DB.Z == i) = tab.fx_L(i);

    rate_gridmin.Z(DB.Z == i) = tab.EvL(i);

    neg_gridmin.Z(DB.Z == i) = tab.vol_err_frac_L(i);

    % K chi and Klp
    % K from chi
    Kchi_gridmin.Z(DB.Z == i) = tab.KL_chi10e6(i);
    
    % Klp
    Klp_gridmin.Z(DB.Z == i) = tab.KlpLchi10e6(i);

    % Morphologic metrics
    
    AreaA.Z(DB.Z == i) = tab.Abkm2(i);

    AspectAsp.Z(DB.Z == i) = tab.asp_ratio(i);

    ReliefR.Z(DB.Z == i) = tab.rel_fac(i);

    SlopeS.Z(DB.Z == i) = tab.slp_fac(i);
    
end

%% PCA analysis

% PCA coloring and indexing by percentile

if strcmp(pca_col,'distance 2 caldera')
    cat_1 = bot_idx';
    cat_3 = top_idx';
elseif strcmp(pca_col,'precipitation')
    cat_1 = top_idx';
    cat_3 = bot_idx';
else strcmp(pca_col,'erosion percentile')
    cat_1 = top_idx_min';
    cat_3 = bot_idx_min';

    cat_1i = top_idx_max';
    cat_3i = bot_idx_max';
end

cat_2 = 1:length(x0);

for i = 1:length(cat_1)
    cat_2(cat_2 == cat_1(i)) = [];
end

for i = 1:length(cat_3)
    cat_2(cat_2 == cat_3(i)) = [];
end

%% Standardize or not for PCA

arr_mat = NaN(size(arr));

pcamin_mat = [];
for i = arr_vars
    
    if strcmp(pca_stand,'Yes')
        mean_arr = mean(arr(:,i));
        max_arr = max(abs(arr(:,i)));
    
        arr_stand = (arr(:,i) - mean_arr)./max_arr;
    
        arr_mat(:,i) = arr_stand;
    elseif strcmp(pca_stand,'No')
        arr_mat(:,i) = arr(:,i);
    else
        disp('Must specify Yes or No to standardization for PCA');
        return
    end
    pcamin_mat = [pcamin_mat,arr_mat(:,i)];
end

[coeff_min,score_min,latent_min,tsquared_min,explained_min] = pca(pcamin_mat);

%% Combine relf_val,slp_val,asp_val valleys in the same quartile into 1
% Vector for each topo metric
val_stats = {relf_val,slp_val,asp_val,[]}; % Same order as arr_vars
comb_stats = {cell(3,1),cell(3,1),cell(3,1),cell(3,1)}; % Same order as arr_vars, combine stats
% Tab all topo stats from each quartile category into 1 vector
for i = 1:length(arr_vars)
    for j = 1:3
        for k = 1:length(eval(['cat_' num2str(j)]))
            if i == 1 || i == 2 % Only slope and area
                comb_stats{1,i}{j,1} = [comb_stats{1,i}{j,1};eval(['val_stats{1,' num2str(i) '}{1,cat_' num2str(j) '(' num2str(k) ')};'])];
            end
        end
    end
end

%% statistical tests on different quartiles using distribution of relstr, slpstr, widstr and area to group quartiles

stats = NaN([length(arr_vars),8]);

for i = 1:length(arr_vars)
    
%%%% HARD CODED FOR QUARTILES ONLY
        q1 = eval(['arr(cat_' num2str(1) ',arr_vars(' num2str(i) '))']);
        q2 = eval(['arr(cat_' num2str(2) ',arr_vars(' num2str(i) '))']);

        [ksh_i,ksp_i] = kstest2(q1,q2);
        [rsp_i,rsh_i] = ranksum(q1,q2);

        stats(i,1) = rsh_i;
        stats(i,2) = rsp_i;
        stats(i,3) = ksh_i;
        stats(i,4) = ksp_i;
        
%%%% HARD CODED FOR QUARTILES ONLY
        q1 = eval(['arr(cat_' num2str(2) ',arr_vars(' num2str(i) '))']);
        q2 = eval(['arr(cat_' num2str(3) ',arr_vars(' num2str(i) '))']);

        [ksh_i,ksp_i] = kstest2(q1,q2);
        [rsp_i,rsh_i] = ranksum(q1,q2);

        stats(i,5) = rsh_i;
        stats(i,6) = rsp_i;
        stats(i,7) = ksh_i;
        stats(i,8) = ksp_i;
        
end

stats_tab = array2table(stats,"VariableNames",...
    ["RS test h (close&mid)","RS test p (close&mid)",...
    "KS test h (close&mid)","KS test p (close&mid)",...
    "RS test h (far&mid)","RS test p (far&mid)",...
    "KS test h (far&mid)","KS test p (far&mid)"]);
stats_tab = [varlab',stats_tab];

%% Regression of climate, volcanism, and erosion

idx_big = tab.Abkm2 > size_thresh;

% Vector of inputs, 2nd input is dependent variable
min_mat = [tab.rain_fac(idx_big)/1,tab.EvL(idx_big),tab.cal_dist(idx_big)/1000];

% Make fit function if you want a polynomial surface
modelfun = @(b,x)b(1) + b(2).*min_mat(:,1) + ...
    b(3).*min_mat(:,3) + b(4).*min_mat(:,1).^2 + ...
    b(5).*min_mat(:,1).*min_mat(:,3) +  b(6).*min_mat(:,3).^2; % Quadratic surface of same function as "fit" above
Xminmat = [min_mat(:,1),min_mat(:,3)];
Yminmat = min_mat(:,2);
fittest = fitnlm(Xminmat,Yminmat,modelfun,[500,0.5,-200,-5e-5,-0.03,20]);

% Make meshgrid for plotting regression
[xq,yq] = meshgrid(linspace(min(min_mat(:,1)),max(min_mat(:,1)),25),...
    linspace(min(min_mat(:,3)),max(min_mat(:,3)),25));
YFIT_min = fittest.Coefficients{1,1} + fittest.Coefficients{2,1}.*xq...
    + fittest.Coefficients{3,1}.*yq + fittest.Coefficients{4,1}.*xq.^2 ...
    + fittest.Coefficients{5,1}.*xq.*yq + fittest.Coefficients{6,1}.*yq.^2; % Quadratic surface

% Calculate point by point values
YFIT_min_pts = fittest.Coefficients{1,1} + fittest.Coefficients{2,1}.*min_mat(:,1)...
    + fittest.Coefficients{3,1}.*min_mat(:,3) + fittest.Coefficients{4,1}.*min_mat(:,1).^2 ...
    + fittest.Coefficients{5,1}.*min_mat(:,1).*min_mat(:,3) + fittest.Coefficients{6,1}.*min_mat(:,3).^2; % Quadratic surface

% Calculate absolute deviation to the point by point values
YFIT_min_pts_Zerr = min_mat(:,2) - YFIT_min_pts;

%% erodibility regressions

cdist_x = min(tab.cal_dist./1000):0.1:max(tab.cal_dist./1000);
precip_x = min(tab.rain_fac)/1000:1/1000:max(tab.rain_fac)/1000;
Ab_x = min(tab.Abkm2):1:max(tab.Abkm2);

fcn = @(b,x) b(1).*exp(b(2).*x) + b(3);
distLmdl = fitnlm(tab.cal_dist./1000,tab.KL_chi10e6,fcn,[1e-6,1e-6,1e-6]);    % K
fL_coeff = [distLmdl.Coefficients{1,1},distLmdl.Coefficients{2,1},distLmdl.Coefficients{3,1}];
fL_ci = coefCI(distLmdl);
fL_mean = fL_coeff(1).*exp(fL_coeff(2).*cdist_x) + fL_coeff(3);
fL_min = fL_ci(1,1).*exp(fL_ci(2,1).*cdist_x) + fL_ci(3,1);
fL_max = fL_ci(1,2).*exp(fL_ci(2,2).*cdist_x) + fL_ci(3,2);
distLmdl_p = coefTest(distLmdl);
[fLrho,fLpval] = corr(tab.cal_dist/1000,tab.KL_chi10e6,'type','Spearman');    % K

[fLrain,gLrain] = fit(tab.rain_fac/1000,tab.KL_chi10e6,'poly1');    % K
fLrain_coeff = coeffvalues(fLrain);
fLrain_ci = confint(fLrain,.95);
fLrain_mean = fLrain_coeff(1).*precip_x + fLrain_coeff(2);
fLrain_min = fLrain_ci(1,1).*precip_x + fLrain_ci(1,2);
fLrain_max = fLrain_ci(2,1).*precip_x + fLrain_ci(2,2);
[fLrainrho,fLrainpval] = corr(tab.rain_fac,tab.KL_chi10e6,'type','Spearman');    % K

rainLmdl = fitlm(tab.rain_fac/1000,tab.KL_chi10e6,'linear');    % K
rainLmdl_p = coefTest(rainLmdl);

s0_x = min(tab.OriginalSlope)/1000:1/1000:max(tab.OriginalSlope)/1000;

[fLs0,gLs0] = fit(tab.OriginalSlope/1000,tab.KL_chi10e6,'poly1');    % K
fLs0_coeff = coeffvalues(fLs0);
fLs0_ci = confint(fLs0,.95);
fLs0_mean = fLs0_coeff(1).*precip_x + fLs0_coeff(2);
fLs0_min = fLs0_ci(1,1).*precip_x + fLs0_ci(1,2);
fLs0_max = fLs0_ci(2,1).*precip_x + fLs0_ci(2,2);
[fLs0rho,fLs0pval] = corr(tab.OriginalSlope,tab.KL_chi10e6,'type','Spearman');    % K

s0Lmdl = fitlm(tab.OriginalSlope/1000,tab.KL_chi10e6,'linear');    % K
s0Lmdl_p = coefTest(s0Lmdl);

%% Plotting parameters
parula9 = flipud(parula(9)); parula9(1,:) = rgb('Khaki');
redbar = autumn(6); redbar(1,:) = rgb('DarkRed'); redbar = flipud(redbar); redbar(1,:) = rgb('Khaki');
redbar5 = autumn(5); redbar5(1,:) = rgb('DarkRed'); redbar5(5,:) = rgb('Khaki');% redbar5 = flipud(redbar5); redbar5(1,:) = rgb('Khaki'); redbar5 = flipud(redbar5);
copper6 = copper(8); copper6(1:2,:) = [];

% Relative point size vector
rel_sz = tab.Abkm2./max(tab.Abkm2)*1000;
rel_sz = rel_sz(idx_big);

%% fprintf text for paper
clc

% Site description

% Range of rainfall atlas annual rainfall
fprintf('Range of annual rainfall = %0.1f - %0.1f\n\n',min(rainatlas.Z(:))/1000,max(rainatlas.Z(:))/1000);

% Range of drainage areas
fprintf('Range of drainage areas = %0.1f - %0.1f\n\n',min(tab.Abkm2),max(tab.Abkm2));

% Results

% best fit theta
fprintf('Best-fit theta = %0.2f\n\n',mn_optimized);

% Percent explained by PCA (only minimum PCA as max not used in paper at
% the moment for PCA)
fprintf('PC1 perc. explained = %0.0f\n',explained_min(1));
fprintf('PC2 perc. explained = %0.0f\n',explained_min(2));
fprintf('PC3 perc. explained = %0.0f\n\n',explained_min(3));

% Multiple linear regression RMSE and R^2
fprintf('Volumetric erosion regression RMSE = %0.0f\n',fittest.RMSE);
fprintf('Volumetric erosion regression p-value = %0.2e\n\n',fittest.ModelFitVsNullModel.Pvalue);

fprintf('Range of volumetric erosion rate = %0.0f - %0.0f\n\n',min(min_mat(:,2)),max(min_mat(:,2)));

% Erodibility regression RMSE and Spearman's rho for caldera distance
fprintf('Erodibility vs dist RMSE = %0.2e\n',distLmdl.RMSE);
fprintf('Erodibility vs dist p = %0.3e\n',distLmdl_p);
fprintf('Erodibility vs dist Spearman Rho = %0.2f\n\n',fLrho);

fprintf('Erodibility vs rain RMSE = %0.2e\n',rainLmdl.RMSE);
fprintf('Erodibility vs rain R2 = %0.2f\n',rainLmdl.Rsquared.Ordinary);
fprintf('Erodibility vs rain p = %0.3e\n',rainLmdl_p);
fprintf('Erodibility vs rain Spearman Rho = %0.2f\n\n',fLrainrho);

% Discussion

rainrange_all = abs(min(tab.rain_fac./1000) - max(tab.rain_fac./1000));
distrange_all = abs(min(tab.cal_dist./1000) - max(tab.cal_dist./1000));

fprintf('Erosion rate decrease of Iao vs Olowalu  = %0.2f\n',(tab.EvL(Olowalu) - tab.EvL(Iao))/tab.EvL(Iao));
fprintf('Erosion rate decrease of Iao vs Ukumehame = %0.2f\n\n',(tab.EvL(Ukumehame) - tab.EvL(Iao))/tab.EvL(Iao));

fprintf('Rain factor range of Iao vs Olowalu compared to full range = %0.2f\n',(tab.rain_fac(Iao) - tab.rain_fac(Olowalu))/rainrange_all/1000);
fprintf('Rain factor range of Iao vs Ukumehame compared to full range = %0.2f\n\n',(tab.rain_fac(Iao) - tab.rain_fac(Ukumehame))/rainrange_all/1000);

fprintf('Honokohau erosion rate is %0.0f perc. that of Iao \n',tab.EvL(Honokohau)/tab.EvL(Iao)*100);
fprintf('Waihee erosion rate is %0.0f perc. that of Iao \n\n',tab.EvL(Waihee)/tab.EvL(Iao)*100);

fprintf('Maximum annual precipitation in Olowalu %0.1f m/yr \n',max(rainatlas.Z(DB.Z == Olowalu)/1000));
fprintf('Maximum annual precipitation in Iao %0.1f m/yr \n',max(rainatlas.Z(DB.Z == Iao)/1000));
fprintf('Maximum annual precipitation in Olowalu/Iao %0.0f perc. \n\n',max(rainatlas.Z(DB.Z == Olowalu)/1000)/max(rainatlas.Z(DB.Z == Iao)/1000)*100);

% ksn_fit K_minfit K_maxfit
fprintf('Ksn chi vs Ksn seg p-value = %0.2e\n',coefTest(ksn_fit));
fprintf('Ksn chi vs Ksn seg RMSE = %0.2f\n',ksn_fit.RMSE);
fprintf('Ksn chi vs Ksn seg R2 = %0.2f\n',ksn_fit.Rsquared.Ordinary);

%% Morphology result maps

morph_col = flipud(hex2rgb(["#993404";"#d95f0e";"#fe9929";"#fed98e";"#ffffd4"]));

figure(1);
tiledlayout(2,2,"TileSpacing","compact");

nexttile
imageschs(dem,ReliefR,'colormap',morph_col,'caxis',[0,1000],'colorbar',true,'ticklabels','none');
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w');
mapshow(basins_MA,'LineStyle','-','FaceColor','none','EdgeColor','k','LineWidth',1);
plot([759384,762384],[2299282,2299282],'k-','LineWidth',3);
set(gca,'XTick',749484);
set(gca,'XTickLabel',['156.6' char(176) 'W']);
set(gca,'YTick',[2302030,2324182]);
set(gca,'YTickLabel',[['20.8' char(176) 'N'];['21.0' char(176) 'N']]);
hold off;
title('Relief');

nexttile
imageschs(dem,SlopeS,'colormap',morph_col,'caxis',[20,45],'colorbar',true,'ticklabels','none');
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w');
mapshow(basins_MA,'LineStyle','-','FaceColor','none','EdgeColor','k','LineWidth',1);
set(gca,'XTick',749484);
set(gca,'YTick',[2302030,2324182]);
hold off;
title('Slope');

nexttile
imageschs(dem,AspectAsp,'colormap',morph_col,'caxis',[0,0.5],'colorbar',true,'ticklabels','none');
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w');
mapshow(basins_MA,'LineStyle','-','FaceColor','none','EdgeColor','k','LineWidth',1);
set(gca,'XTick',749484);
set(gca,'YTick',[2302030,2324182]);
hold off;
title('Aspect ratio');

nexttile
imageschs(dem,AreaA,'colormap',morph_col,'caxis',[0,20],'colorbar',true,'ticklabels','none');
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w');
mapshow(basins_MA,'LineStyle','-','FaceColor','none','EdgeColor','k','LineWidth',1);
set(gca,'XTick',749484);
set(gca,'YTick',[2302030,2324182]);
hold off;
title('Area');

%% PCA figures

PCA_col = [rgb('DarkRed');rgb('DodgerBlue');rgb('Gold')];

figure(2);
if strcmp(pca_col,'valley type')
    tiledlayout(1,2);
elseif strcmp(pca_col,'distance 2 caldera')
    tiledlayout(1,2);
elseif strcmp(pca_col,'precipitation')
    tiledlayout(1,2);
elseif strcmp(pca_col,'erosion percentile')
    tiledlayout(2,2);
end
nexttile;
hold on;
xline(0,'Color',[0.5,0.5,0.5]);yline(0,'Color',[0.5,0.5,0.5]);
scatter(score_min(cat_1,1),score_min(cat_1,2),'m','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(1,:));
scatter(score_min(cat_2,1),score_min(cat_2,2),'b','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(2,:));
scatter(score_min(cat_3,1),score_min(cat_3,2),'c','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(3,:));

for i = 1:length(coeff_min(:,1))

    plot([0,coeff_min(i,1)],[0,coeff_min(i,2)],'k-');
    text(coeff_min(i,1),coeff_min(i,2),varlab{i});

end

xlabel('1st principal component'); ylabel('2nd principal component');
axis square;

nexttile;
hold on;
xline(0,'Color',[0.5,0.5,0.5]);yline(0,'Color',[0.5,0.5,0.5]);
scatter(score_min(cat_1,1),score_min(cat_1,3),'m','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(1,:));
scatter(score_min(cat_2,1),score_min(cat_2,3),'b','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(2,:));
scatter(score_min(cat_3,1),score_min(cat_3,3),'c','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(3,:));

for i = 1:length(coeff_min(:,1))

    plot([0,coeff_min(i,1)],[0,coeff_min(i,3)],'k-');
    text(coeff_min(i,1),coeff_min(i,3),varlab{i});

end

xlabel('1st principal component'); ylabel('3rd principal component');
axis square;

if strcmp(pca_col,'erosion percentile')

nexttile;
hold on;
xline(0,'Color',[0.5,0.5,0.5]);yline(0,'Color',[0.5,0.5,0.5]);
scatter(score_max(cat_1i,1),score_max(cat_1i,2),'m','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(1,:));
scatter(score_max(cat_2i,1),score_max(cat_2i,2),'b','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(2,:));
scatter(score_max(cat_3i,1),score_max(cat_3i,2),'c','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(3,:));

for i = 1:length(coeff_min(:,1))

    plot([0,coeff_min(i,1)],[0,coeff_min(i,2)],'k-');

end

xlabel('1st principal component'); ylabel('2nd principal component');
axis square;

nexttile;
hold on;
xline(0,'Color',[0.5,0.5,0.5]);yline(0,'Color',[0.5,0.5,0.5]);
scatter(score_max(cat_1i,1),score_max(cat_1i,3),'m','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(1,:));
scatter(score_max(cat_2i,1),score_max(cat_2i,3),'b','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(2,:));
scatter(score_max(cat_3i,1),score_max(cat_3i,3),'c','filled','MarkerEdgeColor','none','MarkerFaceColor',PCA_col(3,:));

for i = 1:length(coeff_min(:,1))

    plot([0,coeff_min(i,1)],[0,coeff_min(i,3)],'k-');

end

xlabel('1st principal component'); ylabel('3rd principal component');
axis square;

end

%% PCA distribution figures 

f = figure(3);
tiledlayout(1,length(arr_vars));
% Kernel density cdf
for i = 1:length(arr_vars)

    nexttile;
    hold on;

    if i == 1 || i == 2 % Only slope and relief distributions
        for j = [2,1,3] % Plotting order
            for k = 1:length(eval(['cat_' num2str(j)]))
                eval(['[F1,X1] = ksdensity(val_stats{1,' num2str(i) '}{1,cat_' num2str(j) '(' num2str(k) ')},''Function'',''cdf'');']);
                h = plot(X1,F1,'Color',PCA_col(j,:),'LineWidth',0.5);
                set(h, 'Color', [h.Color, 0.35]);
            end
        end
    end

    % Plot pdfs of mean, min, max, percentile etc topo metrics of all
    % Valleys in a quartile category 
    for j = 1:3
        eval(['[F1,X1] = ksdensity(arr(cat_' num2str(j) ',arr_vars(' num2str(i) ')),''Function'',''cdf'');']);
        plot(X1,F1,'Color',PCA_col(j,:),'LineWidth',2);
    end

    axis square
    xlabel(varlab{i});
    if i == 1
        ylabel('Probability density function');
    end

end
f.Position = [285.4444444444445,517.8888888888889,1116.444444444444,264];

%% Inset erosion rate maps
% Imageschs has trouble plotting correct colors. Plot colors using imagesc
% and use other graphics software to place the transparency on a hillshade

%%%%%%%%% For plotting imagesc when there is an area threshold filter
idx_small = tab.Abkm2 < size_thresh;
idxX_small = find(idx_small);
DB_big = DB;
for i = idxX_small'
    DB_big.Z(find(DB_big.Z == i)) = NaN;
end
DB_big.Z(~isnan(DB_big.Z)) = 1;

parula_long = NaN(255,3);

parula6 = parula(6);

for i = 1:255
    if i == 1
        parula_long(i,:) = [1,1,1];
    elseif i < 42
        parula_long(i,:) = parula6(1,:);
    elseif i < 85
        parula_long(i,:) = parula6(2,:);
    elseif i < 127
        parula_long(i,:) = parula6(3,:);
    elseif i < 170
        parula_long(i,:) = parula6(4,:);
    elseif i < 212
        parula_long(i,:) = parula6(5,:);
    else
        parula_long(i,:) = parula6(6,:);
    end
end

f = figure(4);
tiledlayout(2,1,"TileSpacing","compact");

nexttile;
imageschs(dem,[],'colormap',[0.9,0.9,0.9],'colorbar',true,'ticklabels','none','caxis',[0,1800]);
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w');
title('E_V [t/km2/yr]');

rate_gridminfilt = rate_gridmin; rate_gridminfilt.Z = rate_gridminfilt.Z.*DB_big.Z;

nexttile;
imagesc(rate_gridminfilt);
colormap(parula_long); colorbar;
clim([0,1800]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);

f.Position = [551,361,525,396];

%% Erosion rate regression

f = figure(5);
hold on;
mesh(xq,yq,YFIT_min,'EdgeColor',[0.7,0.7,0.7],'FaceColor','none');
scatter3(min_mat(:,1),min_mat(:,3),min_mat(:,2),rel_sz(idx_big)/2,min_mat(:,2),'filled',...
    'MarkerEdgeColor','k');
xlabel('MAP; m/yr');
ylabel('D_m_i_n; km');
zlabel('E_V; metric ton/km^2/yr');
view(-250,20);
cb = colorbar(); 
ylabel(cb,'E_V; metric ton/km^2/yr','Rotation',270)
cb.Ticks = 0:300:1800;
colormap(parula(6));
clim([0,1800]);
grid on;
box on;
set(gca,'XDir','reverse');
axis square;
xlim([1000,max(min_mat(:,1)) + 500/1000]);
ylim([0,max(min_mat(:,3))]);
zlim([0,1800]);
zticks(0:600:1800);
yticks(0:2:10);
xticklabels(["1","2","3"]);

ytickangle(0);
xtickangle(0);
set(get(gca,'YLabel'),'Rotation',0);
set(get(gca,'XLabel'),'Rotation',0);

f.Position = [289,343,1051,508];

%% Erodibility K vs distance to caldera

winter5 = flipud(winter(6)); winter5(6,:) = rgb('MidnightBlue');
binwid = 3;
bined = 0:binwid:12;
calbin = meanbins(tab.cal_dist(idx_big)/1000,tab.rain_fac(idx_big)/1000,bined);

f = figure(6);

errorbar(tab.cal_dist(idx_big)/1000,tab.KL_chi10e6(idx_big),tab.KL_chi_1s10e6(idx_big),'.','CapSize',0,'Color','k','LineWidth',2);
hold on;
scatter(tab.cal_dist(idx_big)/1000,tab.KL_chi10e6(idx_big),rel_sz/2,tab.rain_fac(idx_big)/1000,'filled');%,'MarkerEdgeColor','k');
plot(cdist_x,fL_mean,'k-');
plot(cdist_x,fL_min,'k--');
plot(cdist_x,fL_max,'k--');
xlabel('Distance to caldera [km]');
ylabel(['K [m^{' num2str(1-m*2) '}/yr]']);
colormap(winter5);
clim([1,4]);

yyaxis right;
errorbar(calbin(:,1),calbin(:,2),calbin(:,2)-calbin(:,3),calbin(:,4)-calbin(:,2),'.','CapSize',10,'Color','b','LineWidth',0.5);
scatter(calbin(:,1),calbin(:,2),30,'w','filled','MarkerEdgeColor','k');
ylim([1,4]);

cb = colorbar();
ylabel(cb,'Mean annual precipitation [m/yr]','Rotation',270);

f.Position =[284,369,1092/2,337];

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';

%% Erodibility K vs precipitation

binwid = 1;
bined = 1:binwid:4;
rainbin = meanbins(tab.rain_fac(idx_big)/1000,tab.cal_dist(idx_big)/1000,bined);

f = figure(7);

errorbar(tab.rain_fac(idx_big)/1000,tab.KL_chi10e6(idx_big),tab.KL_chi_1s10e6(idx_big),'.','CapSize',0,'Color','k','LineWidth',2);
hold on;
scatter(tab.rain_fac(idx_big)/1000,tab.KL_chi10e6(idx_big),rel_sz/2,tab.cal_dist(idx_big)/1000,'filled');%,'MarkerEdgeColor','k');
plot(precip_x,fLrain_mean,'k-');
plot(precip_x,fLrain_min,'k--');
plot(precip_x,fLrain_max,'k--');
xlabel('Precipitation [m/yr]');
ylabel(['K [m^{' num2str(1-m*2) '}/yr]']);
ylim([0,3.5e-6]);
xlim([1,4.5]);
clim([0,10]);

yyaxis right;
errorbar(rainbin(:,1),rainbin(:,2),rainbin(:,2)-rainbin(:,3),rainbin(:,4)-rainbin(:,2),'.','CapSize',10,'Color','r','LineWidth',0.5);
scatter(rainbin(:,1),rainbin(:,2),30,'w','filled','MarkerEdgeColor','k');
set(gca, 'YDir','reverse');
ylim([0,10]);

cb = colorbar();
set(cb,'YDir','reverse' );
colormap(redbar5);

f.Position =[284,369,1092/2,337];

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';

%% Relict slope map

s0col = flipud(hex2rgb(["#980043";"#dd1c77";"#df65b0";"#d7b5d8";"#f1eef6"]));

figure(8);

imageschs(E0,s0,'caxis',[0,20],'colorbar',true,'colormap',s0col,'altitude',70);
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w');
mapshow(basins_MA,'LineStyle','-','FaceColor','none','EdgeColor','k','LineWidth',1);
plot([759384,762384],[2299282,2299282],'k-','LineWidth',3);
hold off;
set(gca,'XTick',749484);
set(gca,'XTickLabel',['156.6' char(176) 'W']);
set(gca,'YTick',[2302030,2324182]);
set(gca,'YTickLabel',[['20.8' char(176) 'N'];['21.0' char(176) 'N']]);

%% Erodibility vs relict slope

f = figure(9);

binwid = 4;
bined = 4:binwid:24;
s0MAPbin = meanbins(tab.OriginalSlope(idx_big),tab.rain_fac(idx_big)/1000,bined);

errorbar(tab.OriginalSlope(idx_big),tab.KL_chi10e6(idx_big),tab.KL_chi_1s10e6(idx_big),'.','CapSize',0,'Color','k','LineWidth',2);
hold on;
scatter(tab.OriginalSlope(idx_big),tab.KL_chi10e6(idx_big),rel_sz/2,tab.rain_fac(idx_big)/1000,'filled');%,'MarkerEdgeColor','k');
xlabel('Relict Slope [deg.]');
ylim([0,3e-6]);
ylabel(['K [m^{' num2str(1-m*2) '}/yr]']);
colormap(winter5);
clim([1,4]);
cb = colorbar();
ylabel(cb,'Mean annual precipitation [m/yr]','Rotation',270);

yyaxis right;
errorbar(s0MAPbin(:,1),s0MAPbin(:,2),s0MAPbin(:,2)-s0MAPbin(:,3),s0MAPbin(:,4)-s0MAPbin(:,2),'.','CapSize',10,'Color','b','LineWidth',0.5)
hold on;
scatter(s0MAPbin(:,1),s0MAPbin(:,2),30,'w','filled','MarkerEdgeColor','k');
ylim([1,4]);

f.Position =[284,369,1092/2,337];
xlim([4,24]);

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
%%
f = figure(10);

binwid = 4;
bined = 4:binwid:24;
s0dminbin = meanbins(tab.OriginalSlope(idx_big),tab.cal_dist(idx_big)/1000,bined);

errorbar(tab.OriginalSlope(idx_big),tab.KL_chi10e6(idx_big),tab.KL_chi_1s10e6(idx_big),'.','CapSize',0,'Color','k','LineWidth',2);
hold on;
scatter(tab.OriginalSlope(idx_big),tab.KL_chi10e6(idx_big),rel_sz/2,tab.cal_dist(idx_big)/1000,'filled');%,'MarkerEdgeColor','k');
xlabel('Relict Slope [deg.]');
ylim([0,3e-6]);
ylabel(['K [m^{' num2str(1-m*2) '}/yr]']);
colormap(redbar5);
clim([0,10]);
cb = colorbar();
set(cb,'YDir','reverse' );
ylabel(cb,'Distance to caldera [km]','Rotation',270);

yyaxis right;
errorbar(s0dminbin(:,1),s0dminbin(:,2),s0dminbin(:,2)-s0dminbin(:,3),s0dminbin(:,4)-s0dminbin(:,2),'.','CapSize',10,'Color','r','LineWidth',0.5)
hold on;
scatter(s0dminbin(:,1),s0dminbin(:,2),30,'w','filled','MarkerEdgeColor','k');
set(gca, 'YDir','reverse');
ylim([0,10]);

xlim([4,24]);
f.Position =[284,369,1092/2,337];

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';

%% K from chi vs Klp

Kcolor = flipud(ttscm('Hawaii'));
Kcolor_spacing = round(linspace(1,255,5));
Kcolor = Kcolor(Kcolor_spacing,:);

f = figure(11);
tiledlayout(1,2,"TileSpacing","tight");
nexttile;
imageschs(dem,Kchi_gridmin,'ticklabels','nice','colormap',Kcolor,'caxis',[0,2.5e-6]);
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w');
mapshow(wai_c_full,'LineStyle',':','FaceColor','none','EdgeColor','w','LineWidth',2);
plot([759384,762384],[2299282,2299282],'k-','LineWidth',3);
patch([759384,760884,762384,760884],[2323182,2327182,2323182,2324182],rgb('black'),'EdgeColor','none');
set(gca,'XTick',749484);
set(gca,'XTickLabel',['156.6' char(176) 'W']);
set(gca,'YTick',[2302030,2324182]);
set(gca,'YTickLabel',[['20.8' char(176) 'N'];['21.0' char(176) 'N']]);
title(['K [m^{' num2str(1-m*2) '}/yr]']);

nexttile;
imageschs(dem,Klp_gridmin,'ticklabels','none','colormap',Kcolor,'caxis',[0,1.25e-6]);
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w');
mapshow(wai_c_full,'LineStyle',':','FaceColor','none','EdgeColor','w','LineWidth',2);
title(['K_l_p [m^{' num2str(1-3*m) '}/yr^{' num2str(1-m) '}]']);
set(gca,'XTick',749484);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[2302030,2324182]);
set(gca,'YTickLabel',[[];[]]);

f.Position = [417.8888888888889,115.6666666666667,843.5555555555554,748.4444444444445/2];

%% compare ksn from chi and segments (main text)

f = figure(12);
plot(ksn_fit);
hold on;
scatter(tab.Cbeta,tab.BA_ksn,50,'filled','k','LineWidth',2);
xlabel(['ksn from \chi [m^' num2str(theta*2,'%0.2f') ']']);
ylabel(['ksn from stream segments [m^' num2str(theta*2,'%0.2f') ']']);
title(['R^2 = ' num2str(ksn_fit.Rsquared.Ordinary,'%0.2f')]);
xlim([80,330]);
ylim([80,330]);
axis square;
hold off;
legend off;
f.Position = [474.7777777777778,474.7777777777778,768.4444444444443,348.0000000000001];

%%

f = figure(13);
tiledlayout(1,2,"TileSpacing","compact");

nexttile;
imageschs(dem,slp,'colormap',redbar,'colorbar',true,'ticklabels','none','caxis',[0,60],'altitude',90);
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w','EdgeAlpha',1,'LineWidth',1);
plot([756652,757652],[2304930,2304930],'k','LineWidth',5);
hold off;
grid off;
set(gca,'XTick',749484);
set(gca,'XTickLabel',['156.6' char(176) 'W']);
set(gca,'YTick',2312946.45);
set(gca,'YTickLabel',['20.9' char(176) 'N']);
xlim([746914.037589336, 757974.927574109]);
ylim([2304473.45256599, 2318132.78865853]);
cb = colorbar(); 
ylabel(cb,'Slope [deg.]','Rotation',270);

nexttile;
imageschs(dem,relf,'colormap',redbar,'colorbar',true,'ticklabels','none','caxis',[0,1200],'altitude',90);
hold on;
mapshow(basins,'LineStyle','-','FaceColor','none','EdgeColor','w','EdgeAlpha',1,'LineWidth',1);
set(gca,'XTick',749484);
set(gca,'YTick',2312946.45);
xlim([746914.037589336, 757974.927574109]);
ylim([2304473.45256599, 2318132.78865853]);
cb = colorbar(); 
ylabel(cb,'1-km radius relief [km]','Rotation',270);

f.Position = [232,246,1210,622];