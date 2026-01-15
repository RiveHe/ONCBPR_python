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
tStart = datetime(2022,7,1);
tEnd   = datetime(2025,7,23);

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
%% STEP 2 – Compute magnitude of completeness (Mc) from magWin

% Clean magnitude vector
magUse = magWin(~isnan(magWin));

if numel(magUse) < 30
    warning('Not enough events in this time window to robustly estimate Mc.');
else
    % Settings for Mc estimation
    dM          = 0.1;   % magnitude bin width
    gofThresh   = 15;    % goodness-of-fit threshold in percent (Wiemer & Wyss use 10–15%)
    Nmin        = 50;    % minimum number of events above Mc for a stable fit

    [Mc, bval, sigma_b, Mc_candidates, R_all] = calc_Mc_Wiemer(magUse, dM, gofThresh, Nmin);

    fprintf('\n===== Magnitude of Completeness Result =====\n');
    fprintf('Time window: %s to %s\n', datestr(tStart), datestr(tEnd));
    fprintf('Number of events in window: %d\n', numel(magUse));
    fprintf('Estimated Mc = %.2f\n', Mc);
    fprintf('b-value (MLE) = %.2f ± %.2f (Shi & Bolt 1982)\n', bval, sigma_b);
    fprintf('Goodness-of-fit threshold = %d%%\n', gofThresh);
    fprintf('===========================================\n\n');

    % Optional: plot magnitude–frequency distribution and Mc
    figure; clf;

    % Histogram of magnitudes (incremental)
    edgesMag   = floor(min(magUse)*10)/10 : dM : ceil(max(magUse)*10)/10;
    centersMag = edgesMag(1:end-1) + dM/2;
    N_inc      = histcounts(magUse, edgesMag);
    N_cum      = fliplr(cumsum(fliplr(N_inc)));   % cumulative N(M >= m)

    % Plot observed cumulative distribution
    semilogy(centersMag, N_cum, 'o', 'MarkerSize', 6, 'DisplayName','Observed cumulative');
    hold on;

    % Overplot theoretical Gutenberg–Richter curve above Mc
    magTh = centersMag(centersMag >= Mc);
    N0    = sum(magUse >= Mc);  % total events above Mc
    N_th  = N0 .* 10.^(-bval .* (magTh - Mc));

    semilogy(magTh, N_th, '-', 'LineWidth', 1.5, 'DisplayName','GR fit (M ≥ Mc)');

    % Mark Mc
    yline(1, '--k', 'HandleVisibility','off'); % just a ref line
    xline(Mc, '--r', sprintf('  Mc = %.2f', Mc), 'LineWidth', 1.2, 'DisplayName','Mc');

    xlabel('Magnitude');
    ylabel('N(M \geq m)');
    title(sprintf('Magnitude–Frequency Distribution (Mc = %.2f, b = %.2f)', Mc, bval));
    legend('Location','southwest');
    grid on;

    % Optional: goodness-of-fit vs Mc candidates
    figure; clf;
    plot(Mc_candidates, R_all, '-o','LineWidth',1.2);
    hold on;
    yline(gofThresh, '--r', sprintf('Threshold = %d%%',gofThresh));
    xlabel('Trial Mc');
    ylabel('Residual (%%)');
    title('Goodness-of-fit vs trial Mc');
    grid on;
end
%%
%% ------------------------------------------------------------------------
%  STEP 0 – user inputs
% -------------------------------------------------------------------------
catFile  = 'H:\Chapter3\LCH\seismic\hypo71.txt';   % full catalog path

% tStart = datetime(2022,7,1);
% tEnd = datetime(2025,9,10);
tStart = datetime(2018,7,1);
tEnd   = datetime(2025,7,1);

% === BPR station coordinates ===
lons = [-129.035593, -129.098678, -129.0819, -129.0988,  -129.1082, -129.1244];
lats = [ 47.958352,  47.948583,  47.9736,   47.9331,     47.924,    47.9599];
names = {'East', 'MEF', 'North', 'South',  'Mothra', 'WF'};

%% STEP 1 – Parse with Depth and MW (Depth comes before MW)
fmt = '%8s %4s %f %f %f %f %f %f %f %*[^\n]';  % ymd, hhmm, sec, latD, latM, lonD, lonM, Depth_km, MW

fid = fopen(catFile,'r');
while ~feof(fid)
    pos = ftell(fid);
    ln  = fgetl(fid);
    if ~isempty(ln) && isstrprop(ln(1),'digit')
        fseek(fid,pos,'bof'); break
    end
end
C = textscan(fid, fmt, 'Delimiter',' ', 'MultipleDelimsAsOne', true);
fclose(fid);

ymd      = C{1};
hhmm     = C{2};
sec      = C{3};
lat_d    = C{4};
lat_m    = C{5};
lon_d    = C{6};
lon_m    = C{7};
Depth_km = C{8};     % depth (km)
mag      = C{9};     % MW

% Convert to decimal degrees
lat = lat_d + lat_m/60;
lon = -(lon_d + lon_m/60);  % West negative

% Build datetime
baseDT = datetime(strcat(ymd, hhmm), 'InputFormat','yyyyMMddHHmm');
dtVec  = baseDT + seconds(sec);

%% ------------------------------------------------------------------------
% NEW STEP – distance from MEF and radial filter (5 km)
% -------------------------------------------------------------------------
% MEF center
mefLon = -129.098678;
mefLat =  47.948583;

R_earth = 6371;  % Earth radius in km

% great-circle distance (haversine)
dlat   = deg2rad(lat - mefLat);
dlon   = deg2rad(lon - mefLon);
latRad = deg2rad(lat);
mefLatRad = deg2rad(mefLat);

a = sin(dlat/2).^2 + cos(mefLatRad).*cos(latRad).*sin(dlon/2).^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));
dist_km = R_earth * c;          % distance from MEF in km

maxR_km = 50;                    % radius = 5 km
inRad   = dist_km <= maxR_km;   % logical index for events within 5 km

% -------------------------------------------------------------------------
% Filter by time window, magnitude (M >= magThres), and radius (<= 5 km)
% -------------------------------------------------------------------------
magThres = 0.6;   % change here if you want a different threshold

inWin   = dtVec >= tStart & dtVec <= tEnd & mag >= magThres & inRad;

latWin  = lat(inWin);
lonWin  = lon(inWin);
magWin  = mag(inWin);
depWin  = Depth_km(inWin);
dtWin   = dtVec(inWin);
distWin = dist_km(inWin);   % distance of kept events (optional)

%% Daily counts (for M ≥ magThres only, within 5 km of MEF)
dayEdges = tStart : caldays(1) : tEnd;          % one-day bins
[N,~]    = histcounts(dtWin, 'BinEdges', dayEdges);
dayMid   = dayEdges(1:end-1)';

xL = [tStart tEnd];   % common x-axis limits

figure('Position',[100 100 1600 600]);
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ------------------------------------------------------------------------
% TOP PANEL – earthquake counts (within 5 km of MEF)
% ------------------------------------------------------------------------
ax1 = nexttile(t,1);
bar(ax1, dayMid, N, 1, 'FaceColor',[0 0 0.5],'EdgeColor','none'); % navy bars
hold(ax1,'on');

% red smoothed curve
winDays   = 7;
N_smooth  = movmean(N, winDays, 'omitnan');
plot(ax1, dayMid, N_smooth, 'r', 'LineWidth', 2);

title(ax1, sprintf('Catalog M \\geq %.1f (within 50 km of MEF)', ...
       magThres));
ylabel(ax1,'Earthquakes per day');
xtickformat(ax1,'yyyy-MM-dd');
grid(ax1,'off'); box(ax1,'on');
ylim(ax1,[0 200]);
xlim(ax1,xL);
ax1.XMinorTick = 'on';
ax1.XAxis.MinorTickValues = tStart : calmonths(1) : tEnd;
ax1.XMinorGrid = 'on';
ax1.FontSize = 18;

% ------------------------------------------------------------------------
% BOTTOM PANEL – time–depth (within 5 km of MEF), marker size = magnitude
% ------------------------------------------------------------------------
ax2 = nexttile(t,2);
if ~isempty(depWin)
    % Scale marker size based on magnitude
    magMin = min(magWin);
    magMax = max(magWin);

    if magMax > magMin
        % Map magnitudes to size range [20, 100] – adjust if you like
        sz = 20 + (magWin - magMin) ./ (magMax - magMin) * 80;
    else
        % All magnitudes are the same
        sz = 50 * ones(size(magWin));
    end

    scatter(ax2, dtWin, depWin, 10, 'filled', 'MarkerFaceAlpha',0.7);
end
set(ax2,'YDir','reverse'); grid(ax2,'on');
xlabel(ax2,'Time'); ylabel(ax2,'Depth (km)');
title(ax2, sprintf('Time–Depth M \\geq %.1f (within 50 km of MEF)', ...
       magThres));
xlim(ax2,xL);
xtickformat(ax2,'yyyy-MM-dd');
ax2.XMinorTick = 'on';
ax2.XAxis.MinorTickValues = tStart : calmonths(1) : tEnd;
ax2.XMinorGrid = 'on';
linkaxes([ax1 ax2],'x');
ax2.FontSize = 18;

%% figure for manuscript

figure('Position',[100 100 1600 600]);
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ------------------------------------------------------------------------
% TOP PANEL – earthquake counts (within 5 km of MEF)
% ------------------------------------------------------------------------
ax1 = nexttile(t,1);
bar(ax1, dayMid, N, 1, 'FaceColor',[0 0 0.5],'EdgeColor','none'); % navy bars
hold(ax1,'on');

% red smoothed curve
winDays   = 7;
N_smooth  = movmean(N, winDays, 'omitnan');
plot(ax1, dayMid, N_smooth, 'r', 'LineWidth', 2);

% title(ax1, sprintf('Catalog M \\geq %.1f (within 20 km of MEF, %s to %s)', ...
%        magThres, datestr(tStart,'yyyy-mm-dd'), datestr(tEnd,'yyyy-mm-dd')));
ylabel(ax1,'Earthquakes per day');
xtickformat(ax1,'yyyy-MM-dd');
grid(ax1,'off'); box(ax1,'on');
ylim(ax1,[0 60]);
xlim(ax1,xL);
ax1.XMinorTick = 'on';
ax1.XAxis.MinorTickValues = tStart : calmonths(1) : tEnd;
ax1.XMinorGrid = 'on';
ax1.FontSize = 18;
% ------------------------------------------------------------------------
% BOTTOM PANEL – time–depth, color = magnitude (optional size = magnitude)
% ------------------------------------------------------------------------
ax2 = nexttile(t,2);

if ~isempty(depWin) && ~isempty(magWin)

    % column vectors
    depWin = depWin(:);
    magWin = magWin(:);
    dtWin  = dtWin(:);

    % remove bad values
    good = ~isnan(depWin) & ~isnan(magWin) & ~isnat(dtWin);
    depWin = depWin(good);
    magWin = magWin(good);
    dtWin  = dtWin(good);

    % ---- OPTIONAL: also scale marker size by magnitude ----
    sMin = 15; sMax = 120;                 % points^2
    magMin = min(magWin);
    magMax = max(magWin);

    if magMax > magMin
        sz = sMin + (magWin - magMin)./(magMax - magMin) * (sMax - sMin);
    else
        sz = (sMin+sMax)/2 * ones(size(magWin));
    end

    % ---- scatter: color by magnitude ----
    sc = scatter(ax2, dtWin, depWin, sz, magWin, 'filled', ...
        'MarkerFaceAlpha',0.8, 'MarkerEdgeColor','none');

    % colormap + colorbar
    colormap(ax2, parula);        % you can change to 'turbo' if you prefer
    cb = colorbar(ax2);
    cb.Label.String = 'Magnitude';
    caxis(ax2, [magMin magMax]);  % lock to data range
end

% Depth direction: shallow at top, deep at bottom
set(ax2,'YDir','reverse'); 
grid(ax2,'on');

xlabel(ax2,'Time');
ylabel(ax2,'Depth (km)');

% title(ax2, sprintf('Time–Depth M \\geq %.1f (within 20 km of MEF, %s to %s)', ...
%        magThres, datestr(tStart,'yyyy-mm-dd'), datestr(tEnd,'yyyy-mm-dd')));

xlim(ax2,xL);
xtickformat(ax2,'yyyy-MM-dd');
ax2.XMinorTick = 'on';
ax2.XAxis.MinorTickValues = tStart : calmonths(1) : tEnd;
ax2.XMinorGrid = 'on';

linkaxes([ax1 ax2],'x');
ax2.FontSize = 18;

