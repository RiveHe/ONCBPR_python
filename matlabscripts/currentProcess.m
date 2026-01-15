%% File list
fps = {
"H:\Chapter3\LCH\ADCP\current\Endeavour_MainEndeavourField_AcousticDopplerCurrentProfiler600kHz_20230701T000000Z_20240101T000000Z-Ensemble3600s_binMapNearest_3beamOn.mat"
"H:\Chapter3\LCH\ADCP\current\Endeavour_MainEndeavourField_AcousticDopplerCurrentProfiler600kHz_20240101T000000Z_20240630T000000Z-Ensemble3600s_binMapNearest_3beamOn.mat"
"H:\Chapter3\LCH\ADCP\current\Endeavour_MainEndeavourField_AcousticDopplerCurrentProfiler600kHz_20240701T000000Z_20241231T000000Z-Ensemble3600s_binMapNearest_3beamOn.mat"
"H:\Chapter3\LCH\ADCP\current\Endeavour_MainEndeavourField_AcousticDopplerCurrentProfiler600kHz_20250101T000000Z_20251015T000000Z-Ensemble3600s_binMapNearest_3beamOn.mat"
"H:\Chapter3\LCH\ADCP\current\Endeavour_MainEndeavourField_AcousticDopplerCurrentProfiler600kHz_20230101T000000Z_20230701T000000Z-Ensemble3600s_binMapNearest_3beamOn.mat"
"H:\Chapter3\LCH\ADCP\current\Endeavour_MainEndeavourField_AcousticDopplerCurrentProfiler600kHz_20220701T000000Z_20221231T000000Z-Ensemble3600s_binMapNearest_3beamOn.mat"};

%% Accumulators
t_all = datetime.empty(1,0);
u_all = []; v_all = []; w_all = [];
pg_all = [];
range_ref = [];

for k = 1:numel(fps)
    S = load(fps{k});
    if isfield(S,'current') && isfield(S.current,'adcp')
        A = S.current.adcp;
    elseif isfield(S,'adcp')
        A = S.adcp;
    else
        error('Could not find *.adcp in file:\n%s', fps{k});
    end

    % Time -> datetime
    t_raw = A.time;
    if isdatetime(t_raw)
        t = t_raw;
    else
        try
            t = datetime(t_raw,'ConvertFrom','datenum');
        catch
            if max(t_raw) > 1e9
                t = datetime(t_raw,'ConvertFrom','posixtime','TimeZone','UTC'); t.TimeZone = '';
            else
                error('Unrecognized time encoding in %s', fps{k});
            end
        end
    end

    % Data
    u = double(A.u);  % (nbin x ntime)
    v = double(A.v);
    if isfield(A,'w'); w = double(A.w); else; w = nan(size(u)); end
    rng = A.range(:);

    if isempty(range_ref)
        range_ref = rng;
    elseif numel(rng) ~= numel(range_ref)
        error('Range bin count differs in file %s', fps{k});
    end

    if isfield(A,'percentGood')
        pg = squeeze(mean(A.percentGood,1)); % nbin x ntime
    else
        pg = nan(size(u));
    end

    % Concatenate along time
    t_all = [t_all, t(:).'];            %#ok<AGROW>
    u_all = [u_all, u];                 %#ok<AGROW>
    v_all = [v_all, v];                 %#ok<AGROW>
    w_all = [w_all, w];                 %#ok<AGROW>
    pg_all = [pg_all, pg];              %#ok<AGROW>
end

% Sort by time & remove duplicates
[ts, idx] = sort(t_all);
u_all = u_all(:,idx); v_all = v_all(:,idx); w_all = w_all(:,idx); pg_all = pg_all(:,idx);
[tsu, ia] = unique(ts);
u_all = u_all(:,ia); v_all = v_all(:,ia); w_all = w_all(:,ia); pg_all = pg_all(:,ia);
ts = tsu;

%% (Optional) enforce exact hourly cadence before Godin
% tt = retime(timetable(ts', (1:size(u_all,2))'), 'hourly','nearest');   % template
% [~, ia2] = ismember(datenum(tt.Time), datenum(ts));
% good = ia2>0;
% if any(good)
%     ts = tt.Time(good);
%     u_all = u_all(:, ia2(good));
%     v_all = v_all(:, ia2(good));
%     w_all = w_all(:, ia2(good));
%     pg_all = pg_all(:, ia2(good));
% end

%% Apply YOUR Godin filter (Z_godin) along time (dim 2)
% Z_godin expects columns = records; it internally transposes when needed.
u_g = Z_godin(u_all);   % returns same size (nbin x ntime), edges become NaN
v_g = Z_godin(v_all);

%% Choose a “seafloor” bin (bottom-mounted up-looking => smallest range)
binIdx = 1;
% if any(isfinite(pg_all(:)))
%     medPG = median(pg_all,2,'omitnan');
%     goodBins = find(medPG >= 20);
%     if ~isempty(goodBins), binIdx = min(goodBins); end
% end
fprintf('Using bin %d (range ≈ %.2f m from transducer)\n', binIdx, range_ref(binIdx));

% Extract filtered near-bed series
u_bot_f = u_g(binIdx,:);
v_bot_f = v_g(binIdx,:);
spd_f   = hypot(u_bot_f, v_bot_f);
time = ts; 
v_bot = v_bot_f;

%% Plot: Godin-filtered near-seafloor currents
figure('Color','w','Position',[80 80 1100 650]);

subplot(3,1,1);
plot(ts, u_bot_f, 'LineWidth',1.2); grid on;
ylabel('u (m/s)');
title(sprintf('Godin-filtered near-seafloor currents'));

subplot(3,1,2);
plot(ts, v_bot_f, 'LineWidth',1.2); grid on;
ylabel('v (m/s)');

subplot(3,1,3);
plot(ts, spd_f, 'LineWidth',1.2); grid on;
ylabel('|U| (m/s)'); xlabel('Time');

%% Optional: Hovmöller of Godin-filtered u,v
figure('Color','w','Position',[80 80 1200 600]);
t_num = datenum(ts);
subplot(1,2,1);
imagesc(t_num, range_ref, u_g); set(gca,'YDir','reverse'); axis tight; colorbar;
title('u_{Godin} (m/s)'); ylabel('range (m)'); datetick('x','keeplimits');

subplot(1,2,2);
imagesc(t_num, range_ref, v_g); set(gca,'YDir','reverse'); axis tight; colorbar;
title('v_{Godin} (m/s)'); ylabel('range (m)'); datetick('x','keeplimits');

%% Optional: Stick plot (filtered) – convert datetime to numeric first
% Inputs expected in workspace:
%   ts        : datetime vector (one per sample)
%   u_bot_f   : Godin-filtered eastward current (m/s)
%   v_bot_f   : Godin-filtered northward current (m/s)

% ---------- settings ----------
t_start = datetime(2022,7,1);
t_end   = datetime(2025,7,1);
ds      = 6;              % plot one vector every 6 hours (change as you like)
ref_speed = 0.05;         % reference arrow = 0.10 m/s (edit to taste)

% ---------- prep ----------
% Keep only data within the window
keep = ts >= t_start & ts <= t_end;
ts_ = ts(keep); u_ = u_bot_f(keep); v_ = v_bot_f(keep);

% Downsample for clarity
idx = 1:ds:numel(ts_);
ts_d = ts_(idx); u_d = u_(idx); v_d = v_(idx);

% Scale by 95th-percentile speed
s = hypot(u_d, v_d);
qscale = max(prctile(s,95), eps);   % avoid divide-by-zero

% X as datenums (quiver doesn't accept datetime)
tx = datenum(ts_d);

% ---------- plot ----------
figure('Color','w','Position',[80 80 1200 300]); hold on; grid on;
quiver(tx, zeros(size(tx)), u_d/qscale, v_d/qscale, 0, ...
       'AutoScale','off', 'MaxHeadSize', 0.4);
yticks([]); ylabel('Current');
title('Stick plot (Godin-filtered; arrows point TOWARD flow)');

% X axis formatting
ax = gca;
ax.XLim = datenum([t_start t_end]);

% nice ticks (every 2 months)
ticks = t_start:calmonths(2):t_end;
ax.XTick = datenum(ticks);
ax.XTickLabel = cellstr(string(ticks, 'yyyy-MM'));

% ---------- reference arrow ----------
% Place a reference vector near the left
x0 = datenum(t_start + days(20));   % 20 days into the window
y0 = 0;                             % baseline
% quiver(x0, y0, ref_speed/qscale, 0, 0, 'Color',[0.2 0.2 0.2], ...
%        'AutoScale','off','MaxHeadSize',0.4, 'LineWidth',1.2);
% text(x0, y0, sprintf('  %0.2f m/s', ref_speed), 'VerticalAlignment','middle');

%(Optional) show north-pointing reference too
quiver(x0+10, y0, 0, ref_speed/qscale, 0, 'Color',[0.2 0.2 0.2], ...
       'AutoScale','off','MaxHeadSize',0.4, 'LineWidth',1.2);
text(x0+10, y0, sprintf('  %0.2f m/s (north)', ref_speed), 'VerticalAlignment','bottom');
%%
%% Plot: all bins — v_g with vertical offsets
figure('Color','w','Position',[100 100 1400 800]); hold on; grid on;

nbin = size(v_g,1);
offset_step = 0.05;                 % m/s vertical spacing between bins (tweak as needed)
C = parula(nbin);                   % colormap across bins (bin 1 -> C(1,:), etc.)

% Choose plotting order: bottom (bin 1) at bottom, upper bins above
bin_order = 1:nbin;                 % or flipud( (1:nbin).' )' to invert

ymin = +inf; ymax = -inf;
for ii = 1:numel(bin_order)
    b = bin_order(ii);
    y = v_g(b,:);                   % ntime
    y = y - mean(y,'omitnan');      % demean per bin (optional, keeps scales comparable)

    yplot = y + (ii-1)*offset_step;
    plot(ts, yplot, 'LineWidth', 1.0, 'Color', C(b,:));

    % baseline for this bin
    yline((ii-1)*offset_step, ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility','off');

    % annotate with range (m) near the left edge where data exists
    ix = find(isfinite(yplot), 1, 'first');
    if ~isempty(ix)
        text(ts(ix) + minutes(30), (ii-1)*offset_step + 0.02, ...
             sprintf('bin %d (%.1f m)', b, range_ref(b)), ...
             'Color', C(b,:), 'FontSize', 9, 'VerticalAlignment','bottom');
    end

    yf = yplot(isfinite(yplot));
    if ~isempty(yf)
        ymin = min(ymin, min(yf));
        ymax = max(ymax, max(yf));
    end
end

xlabel('Time');
ylabel('v (m/s) + offset');
title('Godin-filtered v-velocity at all bins (demeaned, vertically offset)');
set(gca,'FontSize',14,'Box','on');

% Nice limits
if isfinite(ymin) && isfinite(ymax)
    pad = 0.05*(ymax - ymin + eps);
    ylim([ymin - pad, ymax + pad]);
end
% Optional: show a specific time window (uncomment to use)
% xlim([datetime(2023,7,15) datetime(2023,9,23)]);

% Optional colorbar keyed to range
colormap(C);
cb = colorbar; cb.Label.String = 'Bin color (by index)';  % or customize
%%
%% Plot: all bins — u_g with vertical offsets
figure('Color','w','Position',[100 100 1400 800]); hold on; grid on;

nbin = size(u_g,1);
offset_step = 0.05;                 % m/s vertical spacing between bins (tune as needed)
C = parula(nbin);                   % one color per bin
bin_order = 1:nbin;                 % or nbin:-1:1 if you want surface bins on top

ymin = +inf; ymax = -inf;
for ii = 1:numel(bin_order)
    b = bin_order(ii);
    y = u_g(b,:);                   % (ntime)
    y = y - mean(y,'omitnan');      % demean per bin (optional)

    yplot = y + (ii-1)*offset_step;
    plot(ts, yplot, 'LineWidth', 1.0, 'Color', C(b,:));

    % baseline for this bin
    yline((ii-1)*offset_step, ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility','off');

    % annotate with bin & range
    ix = find(isfinite(yplot), 1, 'first');
    if ~isempty(ix)
        text(ts(ix) + minutes(30), (ii-1)*offset_step + 0.02, ...
             sprintf('bin %d (%.1f m)', b, range_ref(b)), ...
             'Color', C(b,:), 'FontSize', 9, 'VerticalAlignment','bottom');
    end

    yf = yplot(isfinite(yplot));
    if ~isempty(yf)
        ymin = min(ymin, min(yf));
        ymax = max(ymax, max(yf));
    end
end

xlabel('Time');
ylabel('u (m/s) + offset');
title('Godin-filtered u-velocity at all bins (demeaned, vertically offset)');
set(gca,'FontSize',14,'Box','on');

if isfinite(ymin) && isfinite(ymax)
    pad = 0.05*(ymax - ymin + eps);
    ylim([ymin - pad, ymax + pad]);
end
% Optional time window:
% xlim([datetime(2023,7,15) datetime(2023,9,23)]);

colormap(C);
cb = colorbar; cb.Label.String = 'Bin color (by index)';

%% Organized 4-panel figure (BPR offsets + u/v/|U|, all time-linked)
f = figure('Color','w','Position',[80 80 1200 800]);
tl = tiledlayout(f,4,1,'TileSpacing','compact','Padding','compact');

% ---------------- Top panel: BPR Godin anomalies with vertical offsets ----------------
ax1 = nexttile; cla(ax1); hold(ax1,'on'); grid(ax1,'on');

% Original lists (North, MEF, South, Mothra, East)
time_list = {time_n,  time_m,  time_s,  time_mo, time_e};
data_list = {northdata_godin, mefdata_godin, southdata_godin, modata_godin, eastdata_godin};
labels    = {'North','MEF','South','Mothra','East'};

% Color mapping by station name (keep station colors consistent)
C2       = lines(numel(labels));
colorFor = @(name) C2(strcmp(labels,name),:);

% East (reference), convert to hPa anomaly once
t_east   = time_e(:);
east_hPa = (eastdata_godin - nanmean(eastdata_godin)) * 100;

% Plot order (top -> bottom) — Mothra, South, MEF, North
plot_order = {'Mothra','South','MEF','North'};
offsets    = [6 4 2 0];   % hPa vertical offsets
order_idx  = cellfun(@(nm) find(strcmp(labels,nm),1), plot_order);  % <— added back

ymin = +inf; ymax = -inf;
for kk = 1:numel(plot_order)
    name = plot_order{kk};
    idx  = order_idx(kk);
    if isempty(idx), continue; end

    t_k = time_list{idx};  if ~isvector(t_k), t_k = t_k(:); end
    x_k = data_list{idx};
    if isempty(t_k) || isempty(x_k), continue; end

    % Convert station to hPa anomaly
    y_k_hPa = (x_k - nanmean(x_k)) * 100;

    % Common time grid with EAST & interpolate both onto it
    t_common = intersect(t_east, t_k);
    if isempty(t_common), continue; end

    e_on = interp1(t_east, east_hPa, t_common, 'linear');
    k_on = interp1(t_k,   y_k_hPa,  t_common, 'linear');

    % Difference (Station − East), demean to show anomaly of the difference
    dP = k_on - e_on;
    dP = dP - nanmean(dP);

    % Plot with vertical offset (keep station color)
    yplot = dP + offsets(kk);
    plot(ax1, t_common, yplot, 'LineWidth', 1.2, 'Color', colorFor(name), ...
        'DisplayName', sprintf('%s - East', name));
    yline(ax1, offsets(kk), ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility','off');

    % Track y-limits
    yf = yplot(isfinite(yplot));
    if ~isempty(yf)
        ymin = min(ymin, min(yf));
        ymax = max(ymax, max(yf));
    end
end

% Axes cosmetics & limits for BPR subplot
ax1.FontSize = 16; ax1.Box = 'on';
if isfinite(ymin) && isfinite(ymax)
    pad = 0.05 * (ymax - ymin + eps);
    ylim(ax1, [ymin - pad, ymax + pad]);
end
legend(ax1,'Location','southoutside','Orientation','horizontal','Box','off');
ylabel(ax1,'\DeltaP (hPa) + offset');
title(ax1,'Godin-Filtered BPR Differences (Station − East)');

% ---------------- Panel 2: u ----------------
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
plot(ax2, ts, u_bot_f, 'LineWidth',1.2);
ylabel(ax2,'u_{Godin} (m/s)');
title(ax2, sprintf('Near-seafloor ADCP (bin %d, range %.1f m)', binIdx, range_ref(binIdx)));
ax2.FontSize = 14; ax2.Box = 'on';

% ---------------- Panel 3: v ----------------
ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on');
plot(ax3, ts, v_bot_f, 'LineWidth',1.2);
ylabel(ax3,'v_{Godin} (m/s)');
ax3.FontSize = 14; ax3.Box = 'on';

% ---------------- Panel 4: |U| ----------------
ax4 = nexttile; hold(ax4,'on'); grid(ax4,'on');
plot(ax4, ts, spd_f, 'LineWidth',1.2);
ylabel(ax4,'|U|_{Godin} (m/s)'); xlabel(ax4,'Time');
ax4.FontSize = 14; ax4.Box = 'on';

% ---------------- Fixed datetime window for ALL panels ----------------
t_start = datetime(2022,8,1);
t_end   = datetime(2025,6,1);
set([ax1,ax2,ax3,ax4],'XLim',[t_start t_end]);

% Link x-axes and avoid duplicate x-tick labels on upper panels
linkaxes([ax1,ax2,ax3,ax4],'x');
xticklabels(ax1,[]); xticklabels(ax2,[]); xticklabels(ax3,[]);
xtickformat(ax4,'yyyy-MM');  % tweak as you like

% Overall figure title (optional)
%title(tl, 'BPR vs. ADCP Currents (Godin-filtered)','FontSize',16,'FontWeight','bold');
%%
% ========= ONE FIGURE: MEF ΔP vs near-bottom v (dual y-axes) =========
f = figure('Color','w','Position',[80 80 1200 500]);  % shorter since it's one panel
ax = axes(f); hold(ax,'on'); grid(ax,'on');

% --- Inputs available from your workspace ---
% time_m, mefdata_godin    % MEF BPR time & series (Godin-filtered)
% time_e, eastdata_godin   % EAST BPR time & series (Godin-filtered)
% ts, v_bot_f              % ADCP time & near-bed v (Godin-filtered)
% range_ref, binIdx        % for the subtitle if you like

% --- Make MEF−EAST ΔP (hPa), demeaned, on a common time base ---
t_east = time_e(:);
t_mef  = time_m(:);
mef_hPa  = (mefdata_godin  - nanmean(mefdata_godin))  * 100;  % anomalies
east_hPa = (eastdata_godin - nanmean(eastdata_godin)) * 100;

t_common = intersect(t_mef, t_east);
mef_on  = interp1(t_mef,  mef_hPa,  t_common, 'linear', NaN);
east_on = interp1(t_east, east_hPa, t_common, 'linear', NaN);

dP_mef = mef_on - east_on;                         % MEF − EAST (hPa)
dP_mef = dP_mef - mean(dP_mef,'omitnan');          % demean

% --- Plot: left y-axis = ΔP (red), right y-axis = v (blue) ---
yyaxis left
p1 = plot(t_common, dP_mef, 'b-', 'LineWidth', 1.4, 'DisplayName','\DeltaP (MEF−East)');
ylabel('\DeltaP (hPa)');
% nice padding
yl = ylim; 
%ylim(yl + [-1 1]*0.05*range(yl));
ylim([-2 4])
yyaxis right
p2 = plot(ts, v_bot_f, 'r-', 'LineWidth', 1.2, 'DisplayName','v_{Godin} (near-bed)');
ylabel('v (m/s)');
ylim([-0.06 0.06])

% --- Cosmetics, legend, titles ---
xlabel('Time');
title(sprintf('MEF − East \\DeltaP vs Near-bed v (bin %d, range %.1f m)', binIdx, range_ref(binIdx)));
%legend([p1 p2], 'Location','southoutside','Orientation','horizontal','Box','off');
set(ax,'FontSize',16,'Box','on');

% --- Fixed datetime window for this figure ---
t_start = datetime(2022,7,1);
t_end   = datetime(2025,7,23);

% t_start = datetime(2023,7,15);
% t_end   = datetime(2023,9,23);
xlim([t_start t_end]);
xtickformat('yyyy-MM');



%% animation 
%% ===== Animation of near-seafloor current + synced ΔP panel =====
% Requires: ts, u_bot_f, v_bot_f  (from your ADCP pipeline)
%           labels, time_list, data_list (for ΔP construction)

% ---- Safety checks ----
varsNeeded = {'ts','u_bot_f','v_bot_f','labels','time_list','data_list'};
for vn = varsNeeded
    if ~exist(vn{1},'var')
        error('Missing variable "%s". Please define it before running this block.', vn{1});
    end
end

% Ensure column vectors for time & components
ts      = ts(:);
u_bot_f = u_bot_f(:);
v_bot_f = v_bot_f(:);

% ---- Output options ----
makeMP4   = true;                       % set true to save MP4
makeGIF   = false;                      % set true to save GIF
mp4_name  = 'MEF_nearbed_current.mp4';
gif_name  = 'MEF_nearbed_current.gif';
fps_out   = 24;                         % video frame rate
max_secs  = 60;                         % ~cap animation length (seconds)
trail_sec = 3*3600;                     % trail length in seconds (e.g., last 3 hours)

% ---- Precompute & sanitize ----
spd_f = hypot(u_bot_f, v_bot_f);
mask  = isfinite(u_bot_f) & isfinite(v_bot_f) & isfinite(ts);
ts    = ts(mask); u_bot_f = u_bot_f(mask); v_bot_f = v_bot_f(mask); spd_f = spd_f(mask);
if isempty(ts)
    error('No finite (u,v) samples to animate.');
end

% Decimate to keep frames reasonable
Nmax     = max_secs * fps_out;
stride   = max(1, ceil(numel(ts) / max(Nmax,1)));
idx_anim = 1:stride:numel(ts);

% Robust axis scale for vector panel
Smax = prctile(spd_f, 99);
if ~isfinite(Smax) || Smax<=0, Smax = max(spd_f); end
L = 1.2 * Smax;

% ---- Figure & layout (2 x 3 tiles) ----
fig = figure('Color','w','Position',[80 80 1200 640]);
tiled = tiledlayout(fig, 2, 3, 'TileSpacing','compact','Padding','compact');

% ===============================
% Panel A: vector (quiver) [2x2]
% ===============================
axV = nexttile(tiled, [2 2]); hold(axV,'on'); box(axV,'on'); grid(axV,'on');
axis(axV,[-L L -L L]); axis(axV,'equal');
xlabel(axV,'East (m/s)'); ylabel(axV,'North (m/s)');
title(axV,'Near-seafloor current vector');

% zero lines
plot(axV, [-L L],[0 0],':','Color',[0.5 0.5 0.5]);
plot(axV, [0 0],[-L L],':','Color',[0.5 0.5 0.5]);

% main arrow + trail (RGB only; no alpha)
hQ = quiver(axV, 0, 0, u_bot_f(idx_anim(1)), v_bot_f(idx_anim(1)), 0, ...
            'LineWidth',2,'MaxHeadSize',0.8,'Color',[0 0 0]);
hTrail = plot(axV, nan, nan, '-', 'LineWidth',1.0, 'Color',[0.6 0.6 0.6]);

% info text (ASCII only)
tInfo = text(axV, 0.98, 0.98, '', 'Units','normalized', ...
             'HorizontalAlignment','right','VerticalAlignment','top', ...
             'BackgroundColor',[1 1 1],'EdgeColor',[0.7 0.7 0.7], ...
             'FontSize',10);

% ===========================================
% Panel B: speed time series (top-right tile)
% ===========================================
axS = nexttile(tiled, 3); hold(axS,'on'); grid(axS,'on'); box(axS,'on');
plot(axS, ts, spd_f, 'LineWidth',1.0);
ylabel(axS,'|U| (m/s)'); title(axS,'Speed');
xlim(axS, [ts(1) ts(end)]);
axes(axS);                                     % for older MATLAB xline behavior
hCursor = xline(ts(idx_anim(1)), '-', 'LineWidth',1.2);

% =========================================
% Panel C: ΔP (Station − East) (bottom-right)
% =========================================
axDP = nexttile(tiled, 6); hold(axDP,'on'); box(axDP,'on'); grid(axDP,'on');

% Build ΔP stacks from provided lists
% Find reference index for 'East' (fallback to 1 if not found)
refIdx = find(strcmp(labels,'East'), 1);
if isempty(refIdx), refIdx = 1; end

nStations = numel(labels);
colors    = lines(nStations);
offset    = 2;     % vertical offset step in hPa
offsetCount = 0;
legendEntries = {};

for sIdx = 1:nStations
    if sIdx == refIdx, continue; end
    t_ref    = time_list{refIdx}(:);
    t_target = time_list{sIdx}(:);
    if isempty(t_ref) || isempty(t_target), continue; end

    % common times (datetime-safe)
    t_common = intersect(t_ref, t_target);
    if isempty(t_common), continue; end

    % series at common times
    d_ref    = interp1(t_ref,    data_list{refIdx}(:), t_common, 'linear', NaN);
    d_target = interp1(t_target, data_list{sIdx}(:),   t_common, 'linear', NaN);

    % difference, demean, convert to hPa, then stack by offset
    d_diff = d_target - d_ref;
    d_diff = d_diff - nanmean(d_diff);
    offsetCount = offsetCount + 1;
    d_plot = d_diff * 100 + (offsetCount-1)*offset;

    plot(axDP, t_common, d_plot, 'LineWidth',1.2, 'Color', colors(sIdx,:));
    legendEntries{end+1} = sprintf('%s - %s', labels{sIdx}, labels{refIdx});
end

title(axDP,'DeltaP (Station - East)');
xlabel(axDP,'Time');
ylabel(axDP,'Delta Pressure (hPa) + offset');
xlim(axDP, [ts(1) ts(end)]);
xtickformat(axDP,'yyyy-MM-dd');

if ~isempty(legendEntries)
    lgdDP = legend(axDP, legendEntries, 'Location','southoutside', ...
                   'Orientation','horizontal','Box','off','Color','none');
end

% synced cursor for ΔP panel
axes(axDP);
hCursorDP = xline(ts(idx_anim(1)), '-', 'LineWidth', 1.2);

% Link x-axes for manual zoom/pan
linkaxes([axS axDP],'x');

% ---- Optional writers ----
if makeMP4
    vw = VideoWriter(mp4_name, 'MPEG-4');
    vw.FrameRate = fps_out;
    vw.Quality   = 95;
    open(vw);
end
if makeGIF
    % no special setup needed
end

% =================
% Animation loop
% =================
for ii = 1:numel(idx_anim)
    k = idx_anim(ii);

    % Update vector panel
    set(hQ, 'UData', u_bot_f(k), 'VData', v_bot_f(k));

    % Update trail (last trail_sec seconds)
    t_now = ts(k);
    keep  = ts >= (t_now - seconds(trail_sec)) & ts <= t_now;
    set(hTrail, 'XData', u_bot_f(keep), 'YData', v_bot_f(keep));

    % Update info & cursors
    dir_deg = mod(atan2d(v_bot_f(k), u_bot_f(k)), 360);
    set(tInfo, 'String', sprintf('%s\nSpeed = %.3f m/s\nDir = %.0f deg (to)', ...
        datestr(t_now, 'yyyy-mm-dd HH:MM:SS'), spd_f(k), dir_deg));
    set(hCursor,   'Value', t_now);
    set(hCursorDP, 'Value', t_now);

    drawnow;

    % Capture frame(s)
    frame = getframe(fig);
    if makeMP4
        writeVideo(vw, frame);
    end
    if makeGIF
        [A,map] = rgb2ind(frame2im(frame), 256);
        if ii==1
            imwrite(A,map,gif_name,'gif','LoopCount',Inf,'DelayTime',1/fps_out);
        else
            imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',1/fps_out);
        end
    end
end

if makeMP4, close(vw); end
disp('Animation complete.');
if makeMP4, fprintf('Saved MP4: %s\n', mp4_name); end
if makeGIF, fprintf('Saved GIF: %s\n', gif_name); end


%%
function smooth = Z_godin(data)
% 24–24–25 Godin filter (Emery & Thomson 1997 Eq. 5.10.37)
% **Use only with hourly data**
% Accepts a vector or a matrix (columns are records). Returns same size.

origSize = size(data);
if origSize(1) < origSize(2), data = data.'; end  % columnize

% Construct 71-point Godin kernel safely (use variable 'f' to avoid name clash)
k = 0:11;
f = NaN(71,1);
f(36:47) = (0.5/(24*24*25))*(1200-(12-k).*(13-k)-(12+k).*(13+k));
k = 12:35;
f(48:71) = (0.5/(24*24*25))*(36-k).*(37-k);
f(1:35)  = flipud(f(37:71));

n   = numel(f);         % 71
pad = floor(n/2);       % 35
N   = size(data,1);
a   = pad + 1;

% If input shorter than kernel, return NaNs of same size
if N < (2*pad + 1)      % < 71 samples
    smooth = NaN(size(data));
    if origSize(1) < origSize(2), smooth = smooth.'; end
    return
end

smooth = zeros(size(data));
for i = 1:size(data,2)
    xi  = data(:,i);
    tmp = conv(xi, f, 'full');           % full conv, then extract
    smooth(:,i) = tmp(a : a+N-1);

    % Safe end padding
    L = min(pad, N);
    smooth(1:L, i) = NaN;
    if N > pad
        smooth(N-pad+1:N, i) = NaN;
    else
        smooth(:, i) = NaN;  % defensive (should be caught above)
    end
end

% Return to original shape if needed
if origSize(1) < origSize(2), smooth = smooth.'; end
end
