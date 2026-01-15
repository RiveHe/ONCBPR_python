%% === Filepaths ===
east_fp  = "C:\Users\River\Downloads\Endeavour_stations\east.mat";
south_fp = "C:\Users\River\Downloads\Endeavour_stations\south.mat";
north_fp = "C:\Users\River\Downloads\Endeavour_stations\north.mat";
mef_fp = "C:\Users\River\Downloads\Endeavour_stations\mef.mat";
mothra_fp = "C:\Users\River\Downloads\Endeavour_stations\mothra.mat";
%westflank_fp = "C:\Users\River\Downloads\Endeavour_stations\westflank.mat";
%mef_southfp = "C:\Users\River\Downloads\Endeavour_stations\mef_south.mat";

% === Load and filter ===
[time_e, eastdata_godin]   = d(east_fp);
[time_s, southdata_godin]  = d(south_fp);
[time_n, northdata_godin]  = d(north_fp);
[time_m, mefdata_godin]  = d(mef_fp);
[time_mo, modata_godin]  = d(mothra_fp);
%[time_w, wfdata_godin]  = load_and_filter(westflank_fp);
%[time_ms, mefsouthdata_godin]  = load_and_filter(mef_southfp); %mef south
%station only have data until 2023-9

%% === Remove Mothra data between 2023-03-19 and 2023-03-23 ===
% === Ensure time_mo is datetime and mask correctly ===
% === Remove Mothra data between 2023-03-18 and 2023-03-23 ===
exclude_start = datetime(2023, 3, 18);
exclude_end   = datetime(2023, 3, 24);
mask = time_mo < exclude_start | time_mo > exclude_end;

time_mo = time_mo(mask);
modata_godin = modata_godin(mask);
%% read HYCOM data
% Define file path
filePath = "H:\Chapter3\Endeavor\HYCOM\HYCOM_predicted_SFP.csv";

% --- Option 1: Read as a table (keeps headers and column names if present)
dataTable = readtable(filePath);

% Display the first few rows
disp(head(dataTable));
% Assumes your table is called dataTable and has:
%   - dataTable.Time (datetime) sampled every 3 hours
%   - dataTable.pressure (numeric)

TT = table2timetable(dataTable,'RowTimes',dataTable.time);

% Apply the classic 24-24-25 hour Godin as three boxcar passes
p1 = movmean(TT.pressure, hours(24), 'SamplePoints', TT.Time, 'Endpoints','fill');
p2 = movmean(p1,            hours(24), 'SamplePoints', TT.Time, 'Endpoints','fill');
p3 = movmean(p2,            hours(25), 'SamplePoints', TT.Time, 'Endpoints','fill');

TT.pressure_godin = p3;

% If you want back a table:
outTable = timetable2table(TT);

%% === Prepare lists ===
time_list = {time_e, time_s, time_n, time_m, time_mo};
data_list = {eastdata_godin, southdata_godin, northdata_godin, mefdata_godin, modata_godin};
labels = {'East', 'South', 'North', 'MEF', 'Mothra'};

% === Set vertical offsets ===
offsets = 0:10:50;

% === Create figure with subplots ===
figure;

% --- Subplot 1: Offset BPR Time Series ---
%subplot(2,1,1);
hold on;

plot(time_e, (eastdata_godin - nanmean(eastdata_godin))*100 + offsets(1), 'LineWidth', 1.2);
plot(time_s, (southdata_godin - nanmean(southdata_godin))*100 + offsets(2), 'LineWidth', 1.2);
plot(time_n, (northdata_godin - nanmean(northdata_godin))*100 + offsets(3), 'LineWidth', 1.2);
plot(time_m, (mefdata_godin - nanmean(mefdata_godin))*100 + offsets(4), 'LineWidth', 1.2);
plot(time_mo, (modata_godin - nanmean(modata_godin))*100 + offsets(5), 'LineWidth', 1.2);
%plot(timedate,bp2_N_nearest + offsets(6), 'LineWidth', 1.2); %caclcualteBP_corpernicus.m
%plot(dataTable.time, (TT.pressure_godin - nanmean(TT.pressure_godin))*100 + offsets(6), 'LineWidth', 1.2);
%plot(time_w, wfdata_godin - nanmean(wfdata_godin) + offsets(6), 'LineWidth', 1.2);

ax = gca; 
ax.FontSize = 16;
ax.Box = 'on';
ylim([-8 60])
legend('East', 'South', 'North', 'MEF', 'Mothra','Copernicus','Location', 'southoutside','Orientation','horizontal','Box','off');
xlabel('Time');
ylabel('Offset Pressure (hPa)');
title('Godin-Filtered and Zero-Mean BPR Time Series (Offset)');
grid on;
%%
figure
hold on

% First line: solid
plot(time_mo, modata_godin - nanmean(modata_godin), 'LineWidth', 1.2,'Color','r');

% Second line: semi-transparent using RGBA (last value is alpha)
c = [0 0.447 0.741]; % default blue color
plot(time_e, eastdata_godin - nanmean(mefdata_godin), ...
    'LineWidth', 1.2, 'Color', [c 0.3]);  % 0.3 = 30% opacity

legend('Mothra Station', 'East Station')
grid on;
%% --- Subplot 2: Pairwise Differences ---
%subplot(2,1,2);
plot_differences_against_reference(time_list, data_list, labels);
%plot_all_station_differences(time_list, data_list, labels);
%%
% === Plot colorful bathymetry map ===
figure;

% Filled contour map (no line styles)
contourf(Lon, Lat, bathy, [-6000:500:0], 'LineStyle', 'none');
colormap(parula);
colorbar;
caxis([-6000 0]);
hold on

% === Overlay contour lines (e.g. -1000, -2000, -3000 m) ===
[C, h] = contour(Lon, Lat, bathy, [-3000 -2000 -1000 0], 'k', 'LineWidth', 1);
clabel(C, h, 'Color', 'k', 'FontSize', 9, 'LabelSpacing', 600);

%% === Plot BPR station markers ===
figure
lons = [-129.035593, -129.098678, -129.0819, -129.0988, -129.0988, -129.1082, -129.1244];
lats = [ 47.958352,  47.948583,  47.9736,  47.9331,   47.9482,   47.924,    47.9599];
names = {'East', 'MEF', 'North', 'South', 'MEF-S', 'Mothra', 'WF'};

plot(lons, lats, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

% Add text labels
for i = 1:length(names)
    text(lons(i)+0.02, lats(i), names{i}, 'FontSize', 9, 'FontWeight', 'bold');
end

% === Plot styling ===
xlabel('Longitude');
ylabel('Latitude');
title('Cascadia Bathymetry with BPR Stations');
axis equal tight;
grid on;
%%
% === Load bathymetry ===
ncfile = "C:\Users\River\Downloads\Endeavour_stations\ncmaps\cascadia.nc";
lon = double(ncread(ncfile, 'x'));  % size [1920]
lat = double(ncread(ncfile, 'y'));  % size [1920]
bathy = double(ncread(ncfile, 'z'));  % size [1920x1920]

% Check orientation and flip if needed
if size(bathy, 1) == length(lon)
    bathy = bathy';  % Make rows = latitude, columns = longitude
end

% === Create meshgrid ===
[Lon, Lat] = meshgrid(lon, lat);  % [lat, lon] shape matches bathy

% === Plot colorful bathymetry ===
figure;
contour_levels = -6000:500:0;
contourf(Lon, Lat, bathy, contour_levels, 'LineStyle', 'none');
colormap(parula);
colorbar;
caxis([-6000 0]);
hold on;

% === Overlay isobaths ===
[C, h] = contour(Lon, Lat, bathy, [0 -500 -1000 -2000 -3000 -4000], 'k');
clabel(C, h, 'Color', 'k', 'FontSize', 8);

% === BPR station coordinates ===
lons = [-129.035593, -129.098678, -129.0819, -129.0988, -129.0988, -129.1082, -129.1244];
lats = [ 47.958352,  47.948583,  47.9736,  47.9331,   47.9482,   47.924,    47.9599];
names = {'East', 'MEF', 'North', 'South', 'MEF-S', 'Mothra', 'WF'};

% === Plot and label BPR stations ===
plot(lons, lats, 'ro', 'MarkerFaceColor', 'r');
for i = 1:length(names)
    text(lons(i)+0.01, lats(i), names{i}, 'FontWeight', 'bold', 'Color', 'w');
end

% === Axis limits and labels ===
xlim([-129.15 -129.00]);
ylim([47.90 48.00]);
xlabel('Longitude');
ylabel('Latitude');
title('Cascadia Bathymetry with BPR Stations');
axis equal;
grid on;
%% plot time series and differences
time_list = {time_e, time_s, time_n, time_m, time_mo};
data_list = {eastdata_godin, southdata_godin, northdata_godin, mefdata_godin, modata_godin};
labels = {'East', 'South', 'North', 'MEF', 'Mothra'};

plot_offset_and_differences_grid(time_list, data_list, labels);
%%
time_list = {time_e, time_s, time_n, time_m, time_mo};
data_list = {eastdata_godin, southdata_godin, northdata_godin, mefdata_godin, modata_godin};
labels = {'East', 'South', 'North', 'MEF', 'Mothra'};

plot_offset_and_east_diff(time_list, data_list, labels);
%% run plotADCP.m first
%% ------------------------------------------------------------------------
%  STEP 0 – user inputs
% -------------------------------------------------------------------------
catFile  = 'H:\Chapter3\LCH\seismic\hypo71.txt';   % full catalog path
%tStart   = datetime(2023,10,01);
%tEnd     = datetime(2023,11,01);
%tStart = datetime(2023,11,23);
%tEnd   = datetime(2024,3,23);
% tStart = datetime(2024,10,25);
% tEnd = datetime(2024,11,25);
% tStart = datetime(2025,01,01);
% tEnd = datetime(2025,01,25);
tStart = datetime(2022,8,1);
tEnd   = datetime(2025,6,1);

% === BPR station coordinates ===
lons = [-129.035593, -129.098678, -129.0819, -129.0988,  -129.1082, -129.1244];
lats = [ 47.958352,  47.948583,  47.9736,  47.9331,   47.924,    47.9599];
names = {'East', 'MEF', 'North', 'South',  'Mothra', 'WF'};
%% STEP 1 – Add magnitude (MW) to format
fmt = '%8s %4s %f %f %f %f %f %f %f %*[^\n]';  % added 1 more %f for MW

fid = fopen(catFile,'r');

while ~feof(fid)
    pos = ftell(fid);
    ln  = fgetl(fid);
    if ~isempty(ln) && isstrprop(ln(1),'digit')
        fseek(fid,pos,'bof');
        break
    end
end

C = textscan(fid, fmt, 'Delimiter',' ', 'MultipleDelimsAsOne',true);
fclose(fid);

ymd    = C{1};
hhmm   = C{2};
sec    = C{3};
lat_d  = C{4};
lat_m  = C{5};
lon_d  = C{6};
lon_m  = C{7};
mag    = C{8};  % magnitude column

% Convert to decimal degrees
lat = lat_d + lat_m / 60;
lon = -(lon_d + lon_m / 60);  % assuming West

% Build datetime
baseDT = datetime(strcat(ymd, hhmm), 'InputFormat','yyyyMMddHHmm');
dtVec  = baseDT + seconds(sec);

% Filter by time window
inWin  = dtVec >= tStart & dtVec <= tEnd;
latWin = lat(inWin);
lonWin = lon(inWin);
magWin = mag(inWin);  % filter magnitudes too
dtWin    = dtVec(inWin);
days = tStart:caldays(1):tEnd;  % or use days(1) if you're using durations
        % daily bin edges
[N,edges] = histcounts(dtWin,days);       % counts per day
dayMid   = edges(1:end-1)';  

%% Plot map with marker size scaled by magnitude
latlim = [47.89, 48];
lonlim = [-129.16, -129];

figure;
worldmap(latlim, lonlim);
% Set font size for latitude and longitude labels
setm(gca, 'FontSize', 15);  % adjust font size as needed


% Normalize magnitude to marker size (adjust scale as needed)
markerSizes = 2 + 5 * (magWin - min(magWin));  % linear scaling

% Plot earthquake epicenters with scaled markers
scatterm(latWin, lonWin, markerSizes, magWin, 'filled');  % marker size ∝ mag^2
colormap('hot');
colorbar;
caxis([min(magWin) max(magWin)]);
title(['Earthquake Epicenters Colored by Magnitude from ', datestr(tStart,'yyyy-mm-dd'), ' to ', datestr(tEnd,'yyyy-mm-dd')], 'FontSize', 18);


% Plot BPR stations
geoshow(lats, lons, 'DisplayType', 'point', ...
        'Marker', 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 10);

% Label BPR stations
for i = 1:length(names)
    textm(lats(i), lons(i) + 0.005, names{i}, ...
          'FontWeight', 'bold', 'Color', 'b', 'FontSize', 14);
end
%%
% Define region and grid size
latlim = [47.89, 48];
lonlim = [-129.16, -129];
dlat = 0.005;  % Grid size in degrees
dlon = 0.005;

% Create grid edges
lat_edges = latlim(1):dlat:latlim(2);
lon_edges = lonlim(1):dlon:lonlim(2);

% 2D histogram: count number of earthquakes in each grid cell
[N, lat_edges, lon_edges] = histcounts2(latWin, lonWin, lat_edges, lon_edges);

% Convert edges to centers for plotting
lat_centers = lat_edges(1:end-1) + dlat/2;
lon_centers = lon_edges(1:end-1) + dlon/2;

% Create meshgrid of lat/lon centers
[LonGrid, LatGrid] = meshgrid(lon_centers, lat_centers);

% Plot map
figure;
worldmap(latlim, lonlim);
setm(gca, 'FontSize', 14)
geoshow(LatGrid, LonGrid, N, 'DisplayType', 'texturemap');
colormap(hot);
colorbar;
title(['Earthquake Counts per Grid Cell from ', ...
    datestr(tStart,'yyyy-mm-dd'), ' to ', datestr(tEnd,'yyyy-mm-dd')], ...
    'FontSize', 16)

% Overlay BPR stations
geoshow(lats, lons, 'DisplayType', 'point', ...
        'Marker', 'o', 'MarkerEdgeColor', 'b', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 8);

for i = 1:length(names)
    textm(lats(i), lons(i) + 0.003, names{i}, ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
end



%% ------------------------------------------------------------------------
%  STEP 3 – figure with 2 rows, 1 column
% -------------------------------------------------------------------------

dayEdges = tStart : caldays(1) : tEnd;          % one-day bins
[N,~]    = histcounts(dtWin, 'BinEdges', dayEdges);
dayMid   = dayEdges(1:end-1)';                  % left edges as x
xL = [datetime(2022,7,1)  datetime(2025,6,1)];   % common x-axis limits

figure('Position',[100 100 1600 600]);           % taller figure
t = tiledlayout(5,1,'TileSpacing','compact','Padding','compact');

% ------------------------------------------------------------------------
% TOP PANEL – earthquake counts
% ------------------------------------------------------------------------
ax1 = nexttile(t,1);
bar(ax1, dayMid, N, 1, 'FaceColor',[0 0 0.5],'EdgeColor','none'); % navy bars
ylim(ax1,[0 2000]);
title(ax1,'Full Catalog','FontWeight','bold');
ylabel(ax1,'Earthquakes per day');
xlim(ax1,xL);

% monthly ticks & YYYY-MM-DD labels
%xticks(ax1, xL(1):calmonths(1):xL(2));
xtickformat(ax1,'yyyy-MM-dd');
%xtickangle(ax1,45);
grid(ax1,'off'); box(ax1,'on');
hold on 


% ------------------------------------------------------------------------
% BOTTOM PANEL – ΔP (each station minus East)
% ------------------------------------------------------------------------
ax2 = nexttile(t,2);  hold(ax2,'on');

refIdx = 1;                       % East = reference
nStations = length(labels);
colors  = lines(nStations);
offset  = 2;
offsetCount = 0;
legendEntries = {};

for sIdx = 1:nStations
    if sIdx == refIdx,  continue;  end

    t_ref    = time_list{refIdx};
    t_target = time_list{sIdx};
    t_common = intersect(t_ref, t_target);
    if isempty(t_common), continue; end

    d_ref    = interp1(t_ref,    data_list{refIdx}, t_common,'linear');
    d_target = interp1(t_target, data_list{sIdx},   t_common,'linear');

    diff = (d_target - d_ref) - nanmean(d_target - d_ref);
    offsetCount = offsetCount + 1;
    diff_plot = diff*100 + (offsetCount-1)*offset;  % scale & offset

    plot(ax2, t_common, diff_plot, 'LineWidth',1.2, 'Color',colors(sIdx,:));
    legendEntries{end+1} = sprintf('%s − East', labels{sIdx});
end

title(ax2,'ΔP (Station − East)');
xlabel(ax2,'Time');
ylabel(ax2,'Δ Pressure (hPa) + offset');
ylim(ax2,[-2 8]);
xlim(ax2,xL);       
xtickformat(ax2,'yyyy-MM-dd');
%xticks(ax2, xL(1):calmonths(1):xL(2));
%xtickangle(ax2,45);
set(ax2,'Box','on');
grid(ax2,'on');

lgd2 = legend(ax2, legendEntries, 'Location','southoutside', ...
              'Orientation','horizontal','Box','off','Color','none');

% Optional: keep the two x-axes locked together when zooming/panning
linkaxes([ax1 ax2],'x');

ax3 = nexttile(t,3);  
% Add offsets for separation
% offset_mef_east    = 0;
% offset_south_east  = 0.03;
% offset_east_mothra = 0.06;
% 
% hold on
% plot(ax3, time, diff_mef_east_demean    + offset_mef_east,    'k-', 'DisplayName','MEF - East')
% plot(ax3, time, diff_south_east_demean  + offset_south_east,  'm-', 'DisplayName','South - East')
% plot(ax3, time, diff_east_mothra_demean + offset_east_mothra, 'b-', 'DisplayName','Mothra - East')
plot(time, mef_data, 'b-')
xlabel('Time')
ylabel('Temperature (°C)')
title('Δ Temperature')
grid on
xlim(ax3,xL) 
% xticks(ax3, xL(1):calmonths(1):xL(2));
xtickformat(ax3,'yyyy-MM-dd');
legend('Location','southoutside', ...
              'Orientation','horizontal','Box','off','Color','none')

ax4 = nexttile(t,4);
plot(time_raw, spd_raw-nanmean(spd_raw), 'b-')
xlabel('Time')
ylabel('Sound speed (m s^{-1})')   % adjust label if this is temperature
title('ADCP Sound-Speed – Main Endeavour Field')
grid on
xlim(ax4,xL) 
%xticks(ax4, xL(1):calmonths(1):xL(2));
xtickformat(ax4,'yyyy-MM-dd');

% ... (your previous code for ax1-ax4 and tiledlayout t)

% ax5 = nexttile(t,5);  % Adds a 5th row to your tiled layout
% 
% plot(ax5, time_ox, ox_raw-nanmean(ox_raw), 'b-')
% xlabel(ax5,'Time')
% ylabel(ax5,'Oxygen (μmol/L)')   % adjust units if needed
% title(ax5,'Oxygen Concentration – East')
% grid(ax5,'on')
% xlim(ax5, xL)
% xtickformat(ax5, 'yyyy-MM-dd');

% ------------------------------------------------------------------------
%  SHADED MASK  •  23-Nov-2023  →  23-Mar-2024
% ------------------------------------------------------------------------
maskStart = datetime(2023,11,23);
maskEnd   = datetime(2024,3,23);
maskCol   = [1 0.6 0.6];      % pink RGB
maskAlpha = 0.25;             % 0 = fully transparent, 1 = opaque

axList = [ax1 ax2 ax3 ax4];   % all axes you want shaded

for ax = axList
    axes(ax);                                   % activate axis
    yl = ylim(ax);                              % current y-range
    p  = patch(ax, [maskStart maskEnd maskEnd maskStart], ...
                    [yl(1)    yl(1)  yl(2) yl(2)], ...
                    maskCol, 'FaceAlpha', maskAlpha, ...
                    'EdgeColor','none');
    % keep the patch in the background so lines/bars stay visible
    uistack(p,'bottom');
end



%%
%% ------------------------------------------------------------------------
%  STEP 4 – correlation analysis  (daily cadence)
% -------------------------------------------------------------------------

% --- make sure both vectors are the same length & shape ---------------
dayMid = dayMid(:);     % force column
N      = N(:);          % force column

if numel(dayMid) ~= numel(N)
    % trim to the shorter of the two – prevents size mismatch
    L = min(numel(dayMid), numel(N));
    dayMid = dayMid(1:L);
    N      = N(1:L);
end

eqTT = timetable(dayMid, N, 'VariableNames', {'EQ'});

% Build a timetable of earthquake counts (already daily).
eqTT = timetable(dayMid, N, 'VariableNames', {'EQ'});

% Pre-allocate result tables
pearson  = NaN(nStations-1,1);
spearman = NaN(nStations-1,1);
ccf_lag  = -60:60;                     % ±60-day window
ccf_mat  = NaN(numel(ccf_lag), nStations-1);

row = 0;
figure('Name','EQ-BPR Scatter','Position',[100 100 1200 300]);
tiledlayout(1,nStations-1,'TileSpacing','compact');

for sIdx = 1:nStations
    if sIdx == refIdx, continue; end           % skip East
    row = row + 1;
    % Step 1: Get common time vector between the target and reference
    t_ref    = time_list{refIdx};
    t_target = time_list{sIdx};
    t_common = intersect(t_ref, t_target);

    if isempty(t_common), continue; end  % skip if no overlap

    % Step 2: Interpolate both time series to common time vector
    d_ref_interp    = interp1(t_ref,    data_list{refIdx}, t_common, 'linear');
    d_target_interp = interp1(t_target, data_list{sIdx},   t_common, 'linear');

    % Step 3: Remove mean and store the result
    diffSeries = (d_target_interp - d_ref_interp) - ...
        nanmean(d_target_interp - d_ref_interp);


    % Put into timetable and DAILY-average to match EQ cadence
    % ------------------------------------------------------------------
    % make sure time and data are column-vectors of identical length
    % ------------------------------------------------------------------
    timeTarget = time_list{sIdx}(:);          % force column
    diffSeries = diffSeries(:);               % force column

    L = min(numel(timeTarget), numel(diffSeries));
    timeTarget = timeTarget(1:L);
    diffSeries = diffSeries(1:L);

    diffTT = timetable(timeTarget, diffSeries, ...
        'VariableNames', {'dP'});

    diffDaily = retime(diffTT, 'daily', 'mean', 'IncludedEdge','left');
    diffDaily.dP = fillmissing(diffDaily.dP,'linear');  % small gaps

    % --- synchronize and de-trend --------------------------------------
    TT = synchronize(eqTT, diffDaily, 'intersection');

    % Optionally remove slow drifts (30-day high-pass)
    % TT.dP = detrend(TT.dP,'linear');     % uncomment if needed
    % TT.EQ = detrend(TT.EQ,'linear');

    % --- (1) zero-lag correlation --------------------------------------
    pearson(row)  = corr(TT.EQ, TT.dP, 'Rows','complete', 'Type','Pearson');
    spearman(row) = corr(TT.EQ, TT.dP, 'Rows','complete', 'Type','Spearman');

    % --- (2) cross-correlation ±60 days --------------------------------
    ccf = xcorr(TT.dP - mean(TT.dP), TT.EQ - mean(TT.EQ), ...
                60, 'coeff');         % dP leads EQ for +lag
    ccf_mat(:,row) = ccf;

    % --- (3) scatter plot ----------------------------------------------
    nexttile; scatter(TT.EQ, TT.dP*100, 10, 'filled');
    lsline;                              % least-squares line
    xlabel('Quakes / day');
    ylabel('\DeltaP (hPa)');
    title(sprintf('%s vs EQ', labels{sIdx}));
end

% ------------------------------------------------------------------------
%  Display lag-0 correlations
% ------------------------------------------------------------------------
T = table(labels(~strcmp(labels,'East'))', pearson, spearman, ...
          'VariableNames',{'Station','Pearson_r','Spearman_\rho'});
disp('Lag-0 correlation (daily means):');
disp(T);

% ------------------------------------------------------------------------
%  Plot cross-correlation functions
% ------------------------------------------------------------------------
figure('Name','Cross-Correlation dP vs EQ','Position',[100 100 700 350]);
hold on
plot(ccf_lag, ccf_mat, 'LineWidth',1.2);
yline(0,'k');
xlabel('Lag (days)   +ve = dP leads EQ');
ylabel('Cross-corr coeff');
legend(labels(~strcmp(labels,'East')), 'Location','eastoutside');
grid on
title('\DeltaP vs Earthquake Count  (daily, ±60-day CCF)');


%%
function plot_offset_and_east_diff(time_list, data_list, labels)
% Plot offset BPR time series (left) and differences from East (right)

nStations = length(labels);
colors = lines(nStations);
y_offset = 10;  % for vertical separation

% === Setup figure with 1 row and 2 columns ===
figure('Position', [100, 200, 1600, 350]);
t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% === Left: Offset Time Series ===
ax1 = nexttile(t, 1); hold(ax1, 'on');

for i = 1:nStations
    data = data_list{i};
    time = time_list{i};
    offset_data = (data - nanmean(data)) * 100 + (i - 1) * y_offset;
    plot(time, offset_data, 'LineWidth', 1.2, 'Color', colors(i, :));
end

lgd = legend(labels, 'Location', 'southoutside');
set(lgd, 'Orientation', 'horizontal', 'Box', 'off', 'Color', 'none');

ax1 = gca;
ax1.FontSize = 16;
ax1.Box = 'on';
xlabel('Time');
ylabel('Offset Pressure (hPa)');
title('Godin-Filtered & Zero-Mean BPR Time Series (Offset)');
ylim([-15 60])
grid on;

% === Right: Differences relative to East ===
ax2 = nexttile(t, 2); hold(ax2, 'on');
refIdx = 1;  % East as reference
legendEntries = {};
offset = 2;  % vertical spacing for visibility
offsetCount = 1;

for targetIdx = 2:nStations  % skip East itself
    if targetIdx == 6
        continue;
    end
    t_ref = time_list{refIdx};
    t_target = time_list{targetIdx};
    t_common = intersect(t_ref, t_target);
    if isempty(t_common), continue; end

    d_ref = interp1(t_ref, data_list{refIdx}, t_common);
    d_target = interp1(t_target, data_list{targetIdx}, t_common);
    diff = d_target - d_ref;
    diff = diff - nanmean(diff);

    offset_diff = diff * 100 + (offsetCount - 1) * offset;
    plot(t_common, offset_diff, 'LineWidth', 1.2, 'Color', colors(targetIdx, :));
    legendEntries{end+1} = sprintf('%s - East', labels{targetIdx});
    offsetCount = offsetCount + 1;
end

ax2 = gca;
ax2.Box = 'on';
ax2.FontSize = 16
xlabel('Time');
ylabel('Δ Pressure + offset');
title('ΔP relative to East');
lgd2 = legend(legendEntries, 'Location', 'southoutside');
set(lgd2, 'Orientation', 'horizontal', 'Box', 'off', 'Color', 'none');
ylim([-2 8])
grid on;
xlim([datetime(2022,7,1) datetime(2025,5,30)])

end

%%
function plot_offset_and_differences_grid(time_list, data_list, labels)
% Plot one offset time series and 5 difference subplots in 3x2 layout

nStations = length(labels);
colors = lines(nStations);
y_offset = 10;  % for vertical separation

% Create figure with 3 columns and 2 rows
figure('Position', [100, 100, 1600, 800]);  % [left, bottom, width, height]
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% === Offset Time Series ===
ax1 = nexttile(t, 1);  % Top-left tile
hold(ax1, 'on');

for i = 1:nStations
    data = data_list{i};
    time = time_list{i};
    offset_data = (data - nanmean(data)) * 100 + (i - 1) * y_offset;
    plot(time, offset_data, 'LineWidth', 1.2, 'Color', colors(i, :));
end

lgd = legend(labels, 'Location', 'southoutside');
set(lgd, 'Orientation', 'horizontal', 'Box', 'off', 'Color', 'none');


xlabel('Time');
ylabel('Offset Pressure (hPa)');
title('Godin-Filtered & Zero-Mean BPR Time Series (Offset)');
ylim([-15 50])
grid on;

% === Pairwise Differences ===
offset = 2.5;
tileIdx = 2;  % Start from tile 2

for refIdx = 1:nStations
    if tileIdx > 6  % Only 6 total tiles
        break;
    end

    ax = nexttile(t, tileIdx); hold(ax, 'on');
    legendEntries = {};
    offsetCount = 1;

    for targetIdx = 1:nStations
        if targetIdx == refIdx, continue; end

        t_ref = time_list{refIdx};
        t_target = time_list{targetIdx};
        t_common = intersect(t_ref, t_target);
        if isempty(t_common), continue; end

        d_ref = interp1(t_ref, data_list{refIdx}, t_common);
        d_target = interp1(t_target, data_list{targetIdx}, t_common);
        diff = d_target - d_ref;
        diff = diff - nanmean(diff);

        offset_diff = diff * 100 + (offsetCount - 1) * offset;
        plot(t_common, offset_diff, 'LineWidth', 1.2, 'Color', colors(offsetCount, :));
        legendEntries{end+1} = sprintf('%s - %s', labels{targetIdx}, labels{refIdx});
        offsetCount = offsetCount + 1;
    end

    title(sprintf('ΔP relative to %s', labels{refIdx}));
    xlabel('Time');
    ylabel('Δ Pressure + offset');
    lgd = legend(legendEntries, 'Location', 'southoutside');
    set(lgd, 'Orientation', 'horizontal', 'Box', 'off', 'Color', 'none');

    %legend(legendEntries, 'Location', 'southeast', 'Box', 'off', 'Color', 'none');
    grid on;
    ylim([-4 10]);

    tileIdx = tileIdx + 1;
end

end


% %%
% function plot_offset_and_differences_grid(time_list, data_list, labels)
% % Plot one offset time series and 5 difference subplots (one per reference station)
% 
% nStations = length(labels);
% colors = lines(nStations);
% y_offset = 10;  % for vertical separation
% 
% % Create figure with wide layout
% figure('Position', [100, 100, 1600, 800]);  % [left, bottom, width, height]
% t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% % === Offset Time Series (wider and shorter) ===
% ax1 = nexttile(t, [2, 1]);  % left panel spans 2 rows only
% hold(ax1, 'on');
% 
% 
% for i = 1:nStations
%     data = data_list{i};
%     time = time_list{i};
%     offset_data = (data - nanmean(data)) * 100 + (i - 1) * y_offset;
%     plot(time, offset_data, 'LineWidth', 1.2, 'Color', colors(i, :));
% end
% 
% lgd = legend(labels, 'Location', 'southeast');  % or any location you used
% set(lgd, 'Box', 'off', 'Color', 'none');
% 
% xlabel('Time');
% ylim([-15 50])
% ylabel('Offset Pressure (hPa)');
% title('Godin-Filtered & Zero-Mean BPR Time Series (Offset)');
% grid on;
% 
% % ===========================================
% % === Pairwise Differences (Each Reference) ==
% % ===========================================
% offset = 2.5;
% colorIdx = 1;
% for refIdx = 1:nStations
%     nexttile; hold on;
%     legendEntries = {};
%     offsetCount = 1;
% 
%     for targetIdx = 1:nStations
%         if targetIdx == refIdx, continue; end
% 
%         % Common time base
%         t_ref = time_list{refIdx};
%         t_target = time_list{targetIdx};
%         t_common = intersect(t_ref, t_target);
%         if isempty(t_common), continue; end
% 
%         d_ref = interp1(t_ref, data_list{refIdx}, t_common);
%         d_target = interp1(t_target, data_list{targetIdx}, t_common);
%         diff = d_target - d_ref;
%         diff = diff - nanmean(diff);
% 
%         offset_diff = diff * 100 + (offsetCount - 1) * offset;
%         plot(t_common, offset_diff, 'LineWidth', 1.2, 'Color', colors(offsetCount, :));
%         legendEntries{end+1} = sprintf('%s - %s', labels{targetIdx}, labels{refIdx});
%         offsetCount = offsetCount + 1;
%     end
% 
%     title(sprintf('ΔP relative to %s', labels{refIdx}));
%     xlabel('Time');
%     ylabel('Δ Pressure + offset');
%     lgd = legend(legendEntries, 'Location', 'southeast');  % or any location you used
%     set(lgd, 'Box', 'off', 'Color', 'none');
%     grid on;
%     ylim([-4 10]);
% end
% 
% end


%% with dedrift step 
function [time_hourly, data_godin_zeroed] = d(filepath)
    % Load .mat file with fields: data.time (datenum), data.dat
    load(filepath);
    time_dt = datetime(data.time, 'ConvertFrom', 'datenum');
    dat = data.dat;

    % ==== CUT DATA BEFORE JULY 2022 ====
    mask = time_dt >= datetime(2022, 8, 1);
    time_dt = time_dt(mask);
    dat = dat(mask);

    % Create timetable and resample to hourly
    TT = timetable(time_dt, dat);
    TT_hourly = retime(TT, 'hourly', 'mean');
    time_hourly = TT_hourly.time_dt;
    data_hourly = TT_hourly.dat;

    % === Remove Linear Drift ===
    valid = ~isnan(data_hourly);
    t_num = datenum(time_hourly);      % convert datetime to datenum for fitting
    p = polyfit(t_num(valid), data_hourly(valid), 1);  % linear fit
    drift = polyval(p, t_num);         % evaluate drift
    data_detrended = data_hourly - drift;

    % === Apply Godin filter ===
    data_godin = Z_godin(data_detrended);

    % === Remove mean (zero-center) ===
    data_godin_zeroed = data_godin - nanmean(data_godin);
end

%% nodedrift step 
% function [time_hourly, data_godin_zeroed] = load_and_filter(filepath)
%     % Load .mat file with fields: data.time (datenum), data.dat
%     load(filepath);
%     time_dt = datetime(data.time, 'ConvertFrom', 'datenum');
%     dat = data.dat;
% 
%     % ==== CUT DATA BEFORE JULY 2022 ====
%     mask = time_dt >= datetime(2022, 7, 1);
%     time_dt = time_dt(mask);
%     dat = dat(mask);
% 
%     % Create timetable and resample to hourly
%     TT = timetable(time_dt, dat);
%     TT_hourly = retime(TT, 'hourly', 'mean');
%     data_hourly = TT_hourly.dat;
% 
%     % Apply Godin filter (must be hourly)
%     data_godin = Z_godin(data_hourly);
% 
%     % Remove mean to zero the time series
%     data_godin_zeroed = data_godin - nanmean(data_godin);
% 
%     time_hourly = TT_hourly.time_dt;
% end
%%
function plot_all_station_differences(time_list, data_list, labels)
% Plot all pairwise station differences on the same figure, with vertical offset
% Inputs:
%   time_list  - cell array of datetime vectors (1 per station)
%   data_list  - cell array of Godin-filtered data vectors (1 per station)
%   labels     - cell array of station names

nStations = length(labels);
nPairs = nchoosek(nStations, 2);
colors = turbo(nPairs);  % more vivid and distinguishable than 'parula'
  % auto color for each pair
colorIdx = 1;
y_offset = 5;         % dbar offset between curves

subplot(212);
hold on;
legendEntries = {};
pairCount = 0;

% Loop through all station pairs
for i = 1:nStations
    for j = i+1:nStations
        % Get common time base
        t1 = time_list{i};
        t2 = time_list{j};
        t_common = intersect(t1, t2);
        if isempty(t_common)
            continue;  % skip if no overlap
        end

        % Interpolate both series to the common time base
        d1 = interp1(t1, data_list{i}, t_common);
        d2 = interp1(t2, data_list{j}, t_common);
        diff = d1 - d2;

        % Remove mean to focus on relative variation
        diff = diff - nanmean(diff);

        % Apply vertical offset
        pairCount = pairCount + 1;
        diff_offset = diff*100 + (pairCount - 1) * y_offset;

        % Plot
        plot(t_common, diff_offset, 'LineWidth', 1.2, 'Color', colors(colorIdx,:));
        legendEntries{end+1} = sprintf('%s - %s', labels{i}, labels{j});
        colorIdx = colorIdx + 1;
    end
end

xlabel('Time');
ylabel('Δ Pressure (hPa) + offset');
title('Godin-Filtered BPR Differences Between Stations (Offset)');
legend(legendEntries, 'Location', 'eastoutside');
grid on;
end
%%
function plot_differences_against_reference(time_list, data_list, labels)
% Plot differences between each station and a reference station in separate figures
% Inputs:
%   time_list  - cell array of datetime vectors (1 per station)
%   data_list  - cell array of Godin-filtered data vectors (1 per station)
%   labels     - cell array of station names

nStations = length(labels);
colors = lines(nStations - 1);  % Enough distinct colors for all comparisons
y_offset = 5;                   % Offset between curves (optional visual spacing)

for refIdx = 1:nStations
    figure('Name', sprintf('Reference: %s', labels{refIdx}), 'Color', 'w');
    hold on;
    legendEntries = {};
    colorIdx = 1;

    for targetIdx = 1:nStations
        if targetIdx == refIdx
            continue;  % skip self
        end

        % Get common time base
        t_ref = time_list{refIdx};
        t_target = time_list{targetIdx};
        t_common = intersect(t_ref, t_target);
        if isempty(t_common)
            continue;
        end

        % Interpolate to common time
        d_ref = interp1(t_ref, data_list{refIdx}, t_common);
        d_target = interp1(t_target, data_list{targetIdx}, t_common);
        diff = d_target - d_ref;

        % Remove mean
        diff = diff - nanmean(diff);

        % Apply vertical offset
        diff_offset = diff * 100 + (colorIdx - 1) * y_offset;

        % Plot
        plot(t_common, diff_offset, 'LineWidth', 1.2, 'Color', colors(colorIdx, :));
        legendEntries{end+1} = sprintf('%s - %s', labels{targetIdx}, labels{refIdx});
        colorIdx = colorIdx + 1;
    end

    title(sprintf('BPR Differences Relative to %s (Godin-filtered)', labels{refIdx}), 'FontWeight', 'bold');
    xlabel('Time');
    ylabel('Δ Pressure (hPa) + offset');
    legend(legendEntries, 'Location', 'eastoutside');
    grid on;
    %xlim([datetime(20)])
end
end
%%
function plot_differences_against_reference_subplot(time_list, data_list, labels)
% Plot differences between each station and a reference station in subplots (same figure)
% Inputs:
%   time_list  - cell array of datetime vectors (1 per station)
%   data_list  - cell array of Godin-filtered data vectors (1 per station)
%   labels     - cell array of station names

nStations = length(labels);
colors = lines(nStations - 1);  % Enough distinct colors
y_offset = 5;

figure('Name', 'All Station Differences by Reference', 'Color', 'w');

for refIdx = 1:nStations
    subplot(nStations, 1, refIdx);
    hold on;
    legendEntries = {};
    colorIdx = 1;

    for targetIdx = 1:nStations
        if targetIdx == refIdx
            continue;
        end

        % Find common time
        t_ref = time_list{refIdx};
        t_target = time_list{targetIdx};
        t_common = intersect(t_ref, t_target);
        if isempty(t_common)
            continue;
        end

        % Interpolate and difference
        d_ref = interp1(t_ref, data_list{refIdx}, t_common);
        d_target = interp1(t_target, data_list{targetIdx}, t_common);
        diff = d_target - d_ref;
        diff = diff - nanmean(diff);
        diff_offset = diff * 100 + (colorIdx - 1) * y_offset;

        plot(t_common, diff_offset, 'LineWidth', 1.2, 'Color', colors(colorIdx, :));
        legendEntries{end+1} = sprintf('%s - %s', labels{targetIdx}, labels{refIdx});
        colorIdx = colorIdx + 1;
    end

    title(sprintf('Reference: %s', labels{refIdx}), 'FontWeight', 'bold');
    ylabel('Δ Pressure (hPa) + offset');
    legend(legendEntries, 'Location', 'eastoutside');
    grid on;
end

xlabel('Time');
end

%%
function [smooth] = Z_godin(data)
% 3/20/2012  Parker MacCready, based on jfilt by Jonathan Lilly 2005
%
% This applies the 24-24-25 Godin Filter to a vector, or a matrix in which
% each COLUMN is a data record.  It differs from the version in
% Z_dasfilt.m in that it makes a single filter shape of length 71, instead
% of using repeated passes of a boxcar. This results in more NaN-padding.
% The other difference is it accepts matrices of column vectors.
% 
% It returns a dataset of the same size you started with, padded whth
% NaN's at the ends.
%
% ** use ONLY with hourly data! **

test = size(data);
% transpose input if necessary
if test(1)<test(2)
    data=data';
end
    
% This is the shape given in Emery and Thomson (1997) Eqn. (5.10.37)
k = [0:11];
filter = NaN * ones(71,1);
filter(36:47) = (0.5/(24*24*25))*(1200-(12-k).*(13-k)-(12+k).*(13+k));
k = [12:35];
filter(48:71) = (0.5/(24*24*25))*(36-k).*(37-k);
filter(1:35) = flipud(filter(37:71));

% % alternatively you can make the shape this way,
% it is identical but less efficient to create
% aa = NaN * ones(25,24,24);
% for ii = [-12:12]
%     aa(ii+13,:,:) = toeplitz([ii:-1:ii-23],[ii:23+ii]);
% end
% X = [-35:35]; % there are 71 values, from -35 to 35, representing indices
% % relative to 0, the offest from the filtered point
% filter = hist(aa(:),X)'; % with distribution N
% filter=filter./sum(filter);

n=length(filter);
smooth=zeros(size(data));
a=round(n/2);
N=size(data,1);
for i=1:size(data,2)
    temp=conv(data(:,i),filter);
    smooth(:,i)=temp(a:a+N-1);
    smooth(1:n,i)=nan*ones(n,1);
    smooth(N-n+1:N,i)=nan*ones(n,1);
end

% return to original shape
if test(1)<test(2)
    smooth=smooth';
end
end
