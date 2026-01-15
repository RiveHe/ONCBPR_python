%% Load ROMS/BPR pair for East & MEF
load("C:\Users\River\Downloads\Endeavour_bottom_pressure_east_mef_2023.mat");
east        = p_b_total_filt(1,:).';
mef         = p_b_total_filt(2,:).';
east_no_atm = p_b_filt(1,:).';
mef_no_atm  = p_b_filt(2,:).';

% Raw time (datenum -> datetime)
time = datetime(date_model_filt, 'ConvertFrom', 'datenum');

%% Load model grid info (used to find south index)
load("C:\Users\River\Downloads\Endeavour_bottom_pressure_2023.mat"); 
south_bpr_lon = -129.0988;
south_bpr_lat = 47.9331;
dist_south = sqrt((lon_model_all - south_bpr_lon).^2 + (lat_model_all - south_bpr_lat).^2);
[~, south_idx] = min(dist_south);

%% Make a clean time vector for the model-aggregated arrays
time_vec = datetime(date_model_filt, 'ConvertFrom', 'datenum');

% --- Ensure p_b_total_filt and p_b_filt have time along rows ---
if size(p_b_total_filt,1) ~= numel(time_vec)
    p_b_total_filt = p_b_total_filt.';   % transpose to [Nt x Nstations]
end
if size(p_b_filt,1) ~= numel(time_vec)
    p_b_filt = p_b_filt.';               % transpose to [Nt x Nstations]
end

% (Optional) unit scaling if your p_b_total_filt was in Pa and you want hPa
% Comment out if not needed.
p_b_total_filt = p_b_total_filt / 10000;
p_b_filt = p_b_filt / 10000; 

% --- Build timetables (now rows must match time_vec) ---
TT          = array2timetable(p_b_total_filt, 'RowTimes', time_vec);
TT_no_atm   = array2timetable(p_b_filt,        'RowTimes', time_vec);

% 2-hourly means
TT_2hour        = retime(TT,        'regular', 'mean', 'TimeStep', hours(2));
TT_2hour_no_atm = retime(TT_no_atm, 'regular', 'mean', 'TimeStep', hours(2));

% Pull arrays back out
time_vec_2hr       = TT_2hour.Time;
p_b_total_filt_2hr = TT_2hour.Variables;        % [N2 x Nstations]
p_b_filt_2hr       = TT_2hour_no_atm.Variables; % [N2 x Nstations]

% Station names (ensure length matches columns)
station_names = {'North','South','North2','Mothra','WF'};
nStations = size(p_b_total_filt_2hr, 2);
if numel(station_names) ~= nStations
    station_names = compose("Sta%02d", 1:nStations);
end

% Align with (east/mef) 1-hour (or original) timestamps
[common_time, idx_2hr, idx_east] = intersect(time_vec_2hr, time);

% Differences relative to East at common times
diff_pb_common  = nan(numel(common_time), nStations);
diff_pb_no_atm  = nan(numel(common_time), nStations);
for i = 1:nStations
    diff_pb_common(:,i) = p_b_total_filt_2hr(idx_2hr,i) - east(idx_east);
    diff_pb_no_atm(:,i) = p_b_filt_2hr(idx_2hr,i)       - east_no_atm(idx_east);  % <-- fixed typo
end

% MEF - East (with and without atm)
diff_mef        = mef(idx_east)        - east(idx_east);
diff_mef_no_atm = mef_no_atm(idx_east) - east_no_atm(idx_east);

%% Plot
figure
subplot(2,1,1); hold on
cols = lines(nStations);
for i = 1:nStations
    if (i == south_idx && south_idx <= nStations) || i == 3
        continue  % Skip south station (if within range) and i==3 as in your code
    end
    plot(common_time, (diff_pb_common(:,i) - nanmean(diff_pb_common(:,i)))*100, ...
        'LineWidth',1.2, 'Color',cols(i,:), 'DisplayName',[station_names{i} ' - East']);
end
plot(common_time, (diff_mef - nanmean(diff_mef))*100, 'k-', 'LineWidth',1.2, 'DisplayName','MEF - East')
xlim([datetime(2023,07,15), datetime(2023,09,23)])
xlabel('Time'); ylabel('BP Difference (hPa)')
title('Bottom Pressure Difference Relative to East Station (ROMS with atm)')
legend('Location','best'); grid on

subplot(2,1,2); hold on
for i = 1:nStations
    if (i == south_idx && south_idx <= nStations) || i == 3
        continue
    end
    plot(common_time, (diff_pb_no_atm(:,i) - nanmean(diff_pb_no_atm(:,i)))*100, ...
        'LineWidth',1.2, 'Color',cols(i,:), 'DisplayName',[station_names{i} ' - East']);
end
plot(common_time, (diff_mef_no_atm - nanmean(diff_mef_no_atm))*100, 'k-', 'LineWidth',1.2, 'DisplayName','MEF - East')
xlim([datetime(2023,07,15), datetime(2023,09,23)])
xlabel('Time'); ylabel('BP Difference (hPa)')
title('Bottom Pressure Difference Relative to East Station (ROMS no atm)')
legend('Location','best'); grid on

%%
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


%% plot BPR data

% Interpolate other stations to match south station's time grid
eastdata_on_south = interp1(time_e, eastdata_godin, time_s, 'linear', NaN);
northdata_on_south = interp1(time_n, northdata_godin, time_s, 'linear', NaN);
mefdata_on_south = interp1(time_m, mefdata_godin, time_s, 'linear', NaN);
modata_on_south = interp1(time_mo, modata_godin, time_s, 'linear', NaN); % already masked

% Optionally, check alignment:
% figure; plot(time_s, eastdata_on_south - southdata_godin)
%%
diff_east_south  = southdata_godin  - eastdata_on_south;
diff_north_south = northdata_on_south - eastdata_on_south;
diff_mef_south   = mefdata_on_south   - eastdata_on_south;
diff_mothra_south = modata_on_south   - eastdata_on_south;

%%
figure
hold on

offset_step = 1; % vertical offset in hPa
pair_names = {'South - East','North - East','MEF - East','Mothra - East'};
bpr_data = {diff_east_south, diff_north_south, diff_mef_south, diff_mothra_south};
roms_data = {diff_pb_common(:,2), diff_pb_common(:,1), diff_mef, diff_pb_common(:,4)};
colors = {'r','g','b','m'};  % colors for BPR lines
roms_colors = {'r--','g--','b--','m--'}; % dashed for ROMS

for i = 1:length(pair_names)
    offset = (i-1)*offset_step;

    % BPR line
    plot(time_s, (bpr_data{i}-nanmean(bpr_data{i}))*100 + offset, ...
        colors{i}, 'LineWidth',1.5,'DisplayName',[pair_names{i} ' BPR']);

    % BPR ROMS line (dashed)
    plot(common_time, (roms_data{i}-nanmean(roms_data{i}))*100 + offset, ...
        roms_colors{i}, 'LineWidth',1.5,'DisplayName',[pair_names{i} ' ROMS']);
end

ax = gca;
ax.FontSize = 17;
xlim([datetime(2023,07,15), datetime(2023,09,23)])
xlabel('Time')
ylabel('Godin-filtered BP Difference (hPa)')
title('Bottom Pressure Difference Relative to East Station (BPR & ROMS)')
grid on
legend('Location','eastoutside')
hold off
%% compare BPR and ROMS (differences and demeaned)
% === Two-panel figure: (1) Differences vs East (your original idea)
%                        (2) Absolute station series (BPR vs ROMS) ===

figure('Position',[100 100 1600 700])
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% -------------------------
% TOP: Differences vs East
% -------------------------
nexttile; hold on

offset_step = 1; % vertical offset in hPa
pair_names = {'South - East','North - East','MEF - East','Mothra - East'};
bpr_data   = {diff_east_south,  diff_north_south,  diff_mef_south,  diff_mothra_south};
roms_data  = {diff_pb_common(:,2), diff_pb_common(:,1), diff_mef, diff_pb_common(:,4)};
colors     = {'r','g','b','m'};      % BPR solid
roms_styles= {'r--','g--','b--','m--'}; % ROMS dashed

for i = 1:numel(pair_names)
    off = (i-1)*offset_step;

    % BPR
    plot(time_s, (bpr_data{i}-nanmean(bpr_data{i}))*100 + off, ...
        colors{i}, 'LineWidth',1.5, 'DisplayName',[pair_names{i} ' BPR']);

    % ROMS
    plot(common_time, (roms_data{i}-nanmean(roms_data{i}))*100 + off, ...
        roms_styles{i}, 'LineWidth',1.5, 'DisplayName',[pair_names{i} ' ROMS']);
end

ax = gca; ax.FontSize = 16;
xlim([datetime(2023,07,15) datetime(2023,09,23)])
xlabel('Time')
ylabel('Godin-filtered BP Difference (hPa)')
title('Bottom Pressure Difference Relative to East Station (BPR & ROMS)')
grid on
legend('Location','eastoutside')
hold off

% ---------------------------------------------------------
% BOTTOM: Absolute station series (BPR vs ROMS), not diffs
% ---------------------------------------------------------
% Expect the following (Godin-filtered) absolute series to exist:
%   east_bpr, south_bpr, north_bpr, mef_bpr, mothra_bpr
%   east_roms, south_roms, north_roms, mef_roms, mothra_roms
% Each is a column vector aligned with time_s (BPR) or common_time (ROMS)

% Plot BPR & ROMS together per station; each next station sits ABOVE the previous
nexttile; hold on

station_names = {'East','South','North','MEF','Mothra'};
bpr_abs  = {eastdata_godin, southdata_godin, northdata_godin, mefdata_godin, modata_godin};
time_list = {time_e,        time_s,          time_n,          time_m,        time_mo};
roms_abs = {east(idx_east), p_b_total_filt_2hr(idx_2hr,2),    p_b_total_filt_2hr(idx_2hr,1), ...
            mef(idx_east),  p_b_total_filt_2hr(idx_2hr,4)};

abs_colors     = {'k','r','g','b','m'};           % solid for BPR
abs_roms_style = {'k--','r--','g--','b--','m--'}; % dashed for ROMS

offset_step = 8.0;  % larger vertical spacing so stations do NOT overlap

for i = 1:numel(station_names)
    off = (i-1)*offset_step;  % each station at a different vertical baseline

    % BPR and ROMS at the same vertical offset for the same station
    bpr_plot  = (bpr_abs{i}  - nanmean(bpr_abs{i}))*100 + off;
    roms_plot = (roms_abs{i} - nanmean(roms_abs{i}))*100 + off;

    plot(time_list{i}, bpr_plot,  abs_colors{i},     'LineWidth',1.5, ...
        'DisplayName',[station_names{i} ' BPR']);
    plot(common_time, roms_plot, abs_roms_style{i}, 'LineWidth',1.5, ...
        'DisplayName',[station_names{i} ' ROMS']);
end

ax2 = gca; ax2.FontSize = 16;
xlim([datetime(2023,07,15) datetime(2023,09,23)])
xlabel('Time')
ylabel('Godin-filtered BP (hPa)')
title('Demeaned Bottom Pressure by Station (BPR & ROMS)')
grid on
legend('Location','eastoutside')
hold off

%% plot three panels for endeavour proposal figure D
addpath("C:\Users\River\Downloads\pmtmPH_cmtm")
% === Three-panel figure: (1) Differences vs East
%                         (2) Absolute station series (BPR vs ROMS)
%                         (3) Spectra: ROMS vs BPR at each station (offset)
figure('Position',[100 100 1600 850])
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% -------------------------
% LEFT: Differences vs East
% -------------------------
nexttile; hold on
offset_step = 1; % vertical offset in hPa
pair_names = {'South - East','North - East','MEF - East','Mothra - East'};
bpr_data   = {diff_east_south,  diff_north_south,  diff_mef_south,  diff_mothra_south};
roms_data  = {diff_pb_common(:,2), diff_pb_common(:,1), diff_mef, diff_pb_common(:,4)};
colors     = {'r','g','b','m'};      % BPR solid
roms_styles= {'r--','g--','b--','m--'}; % ROMS dashed

for i = 1:numel(pair_names)
    off = (i-1)*offset_step;
    plot(time_s, (bpr_data{i}-nanmean(bpr_data{i}))*100 + off, ...
        colors{i}, 'LineWidth',1.5, 'DisplayName',[pair_names{i} ' BPR']);
    plot(common_time, (roms_data{i}-nanmean(roms_data{i}))*100 + off, ...
        roms_styles{i}, 'LineWidth',1.5, 'DisplayName',[pair_names{i} ' ROMS']);
end
ax = gca; ax.FontSize = 16;
xlim([datetime(2023,07,15) datetime(2023,09,23)])
xlabel('Time'); ylabel('Godin-filtered BP Difference (hPa)')
title(sprintf('Bottom Pressure Difference Relative to East Station \nWith Atmospheric Pressure'))
grid on
lg1 = legend('NumColumns',2,'Location','southoutside');  % two-row legend at bottom
hold off

% ---------------------------------------------------------
% MIDDLE: Absolute station series (BPR vs ROMS), not diffs
% ---------------------------------------------------------
nexttile; hold on
station_names = {'East','South','North','MEF','Mothra'};
bpr_abs  = {eastdata_godin, southdata_godin, northdata_godin, mefdata_godin, modata_godin};
time_list = {time_e,        time_s,          time_n,          time_m,        time_mo};
roms_abs = {east(idx_east), p_b_total_filt_2hr(idx_2hr,2),    p_b_total_filt_2hr(idx_2hr,1), ...
            mef(idx_east),  p_b_total_filt_2hr(idx_2hr,4)};
abs_colors     = {'k','r','g','b','m'};           % solid for BPR
abs_roms_style = {'k--','r--','g--','b--','m--'}; % dashed for ROMS
offset_step_abs = 8.0;  % vertical spacing (hPa)

for i = 1:numel(station_names)
    off = (i-1)*offset_step_abs;
    bpr_plot  = (bpr_abs{i}  - nanmean(bpr_abs{i}))*100 + off;
    roms_plot = (roms_abs{i} - nanmean(roms_abs{i}))*100 + off;
    plot(time_list{i}, bpr_plot,  abs_colors{i},     'LineWidth',1.5, ...
        'DisplayName',[station_names{i} ' BPR']);
    plot(common_time, roms_plot, abs_roms_style{i}, 'LineWidth',1.5, ...
        'DisplayName',[station_names{i} ' ROMS']);
end
ax2 = gca; ax2.FontSize = 16;
xlim([datetime(2023,07,15) datetime(2023,09,23)])
xlabel('Time'); ylabel('Godin-filtered BP (hPa)')
title(sprintf('Demeaned Bottom Pressure by Station \nWith Atmospheric Pressure'))
grid on
lg2 = legend('NumColumns',2,'Location','southoutside');
hold off

% ----------------------------------------------------------------
% ----------------------------------------------------------------
% RIGHT: Spectral comparison (ROMS vs BPR) with vertical offsets
%       -> now using multitaper pmtmPH_ori (uniform sampling)
% ----------------------------------------------------------------
nexttile; hold on
offset_step_spec = 1.0;     % vertical spacing between stations (log10 power units)
t_start = datetime(2023,07,15);
t_end   = datetime(2023,09,23);
nw = 1;                     % time-bandwidth product for MTM

for i = 1:numel(station_names)
    off = (i-1)*offset_step_spec;

    % ===== BPR (hourly) =====
    maskTimeB = time_list{i} >= t_start & time_list{i} <= t_end;
    xB = bpr_abs{i}(maskTimeB);

    % Fill small gaps to keep MTM stable; demean
    xB = fillmissing(xB,'linear','EndValues','nearest');
    xB = xB - nanmean(xB);
    nB = numel(xB);

    if nB >= 8
        % dt in DAYS as a scalar double
        dt_bpr_days = days(median(diff(time_list{i}(maskTimeB))));
        if ~isfinite(dt_bpr_days) || dt_bpr_days<=0, dt_bpr_days = 1/24; end

        % MTM spectrum (pmtmPH_ori returns [P, f] where f is cycles/day if dt is days)
        [PB, fB] = pmtmPH_ori(xB, dt_bpr_days, nw, 0, nB);
        plot(fB, log10(PB)+off, abs_colors{i}, 'LineWidth',1, ...
            'DisplayName',[station_names{i} ' BPR']);
    end

    % ===== ROMS (2-hourly) =====
    maskTimeR = common_time >= t_start & common_time <= t_end;
    xR = roms_abs{i}(maskTimeR);

    xR = fillmissing(xR,'linear','EndValues','nearest');
    xR = xR - nanmean(xR);
    nR = numel(xR);

    if nR >= 8
        dt_roms_days = days(hours(2));  % 2 hours in days, as scalar double
        [PR, fR] = pmtmPH_ori(xR, dt_roms_days, nw, 0, nR);
        plot(fR, log10(PR)+off, abs_roms_style{i}, 'LineWidth',1, ...
            'DisplayName',[station_names{i} ' ROMS']);
    end
end

ax3 = gca; ax3.FontSize = 16;
xlim([0 1])                      % frequency range in cycles/day
xlabel('Frequency (cycles/day)')
ylabel('log_{10} Power (offset)')
title(sprintf('Spectral Comparison by Station \nWith Atmospheric Pressure'))
grid on
lg3 = legend('NumColumns',2,'Location','southoutside');
hold off

%% no atmospheric pressure
addpath("C:\Users\River\Downloads\pmtmPH_cmtm")

% === Define common colormap ===
nStations = 5;
colors = lines(nStations);
abs_colors = colors(1:nStations, :);  % Solid lines (BPR)
abs_lineStyles = repmat({'-'}, 1, nStations);  % For BPR
roms_lineStyles = repmat({'--'}, 1, nStations);  % For ROMS

% === Station names and data ===
station_names = {'East','South','North','MEF','Mothra'};
bpr_abs  = {eastdata_godin, southdata_godin, northdata_godin, mefdata_godin, modata_godin};
time_list = {time_e,        time_s,          time_n,          time_m,        time_mo};
roms_abs = {east_no_atm(idx_east), p_b_filt_2hr(idx_2hr,2),    p_b_filt_2hr(idx_2hr,1), ...
            mef_no_atm(idx_east),  p_b_filt_2hr(idx_2hr,4)};

% === Difference data ===
pair_names = {'South - East','North - East','MEF - East','Mothra - East'};
bpr_data   = {diff_east_south,  diff_north_south,  diff_mef_south,  diff_mothra_south};
roms_data  = {diff_pb_no_atm(:,2), diff_pb_no_atm(:,1), diff_mef_no_atm, diff_pb_no_atm(:,4)};

% === Plot setup ===
figure('Position',[100 100 1600 950])
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% --------------------------------------
% LEFT: Differences relative to East
% --------------------------------------
% --------------------------------------
% LEFT: Differences relative to East
% --------------------------------------
nexttile; hold on
offset_step = 1; % vertical offset in hPa

for i = 1:4
    off = (i-1)*offset_step;
    bpr_interp = interp1(time_s, bpr_data{i}, common_time, 'linear', 'extrap');
    bpr_diff = (bpr_interp - nanmean(bpr_interp)) * 100;
    roms_diff = (roms_data{i}-nanmean(roms_data{i}))*100;
    residual = roms_diff - bpr_diff;

    % === BPR difference ===
    plot(common_time, bpr_diff + off, ...
        'Color', abs_colors(i+1, :), 'LineStyle', '-', 'LineWidth',1.5, ...
        'DisplayName',[pair_names{i} ' BPR']);

    % === ROMS difference ===
    plot(common_time, roms_diff + off, ...
        'Color', abs_colors(i+1, :), 'LineStyle', '--', 'LineWidth',1.5, ...
        'DisplayName',[pair_names{i} ' ROMS']);

    % === Residual ===
    h_res = plot(common_time, residual + off, ...
        'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth',1.5);
    
    if i == 1
        h_res.DisplayName = 'Residual (ROMS - BPR)';
    else
        h_res.Annotation.LegendInformation.IconDisplayStyle = 'off';  % Exclude from legend
    end
end



ax = gca; ax.FontSize = 16;
xlim([datetime(2023,07,15) datetime(2023,09,23)])
xlabel('Time'); ylabel('Godin-filtered BP Difference (hPa)')
title('BP Difference Relative to East (BPR, ROMS, Residual)')
legend('NumColumns',2,'Location','southoutside');
grid on

% ---------------------------------------------------------
% MIDDLE: Absolute BPR vs ROMS with offsets
% ---------------------------------------------------------
nexttile; hold on
offset_step_abs = 8.0;

for i = 1:nStations
    off = (i-1)*offset_step_abs;
    bpr_plot  = (bpr_abs{i}  - nanmean(bpr_abs{i}))*100 + off;
    roms_plot = (roms_abs{i} - nanmean(roms_abs{i}))*100 + off;

    plot(time_list{i}, bpr_plot, 'Color', abs_colors(i,:), ...
        'LineStyle', '-', 'LineWidth',1.5, ...
        'DisplayName',[station_names{i} ' BPR']);
    plot(common_time, roms_plot, 'Color', abs_colors(i,:), ...
        'LineStyle', '--', 'LineWidth',1.5, ...
        'DisplayName',[station_names{i} ' ROMS']);
end

ax2 = gca; ax2.FontSize = 16;
xlim([datetime(2023,07,15) datetime(2023,09,23)])
xlabel('Time'); ylabel('Godin-filtered BP (hPa)')
title('Demeaned Bottom Pressure by Station')
legend('NumColumns',2,'Location','southoutside');
grid on

% ---------------------------------------------------------
% RIGHT: Spectral comparison (ROMS vs BPR) with offsets
% ---------------------------------------------------------
nexttile; hold on
offset_step_spec = 3.0;
t_start = datetime(2023,07,15);
t_end   = datetime(2023,09,23);
nw = 1;  % time-bandwidth product

for i = 1:nStations
    off = (i-1)*offset_step_spec;

    % === BPR spectrum ===
    maskTimeB = time_list{i} >= t_start & time_list{i} <= t_end;
    xB = fillmissing(bpr_abs{i}(maskTimeB), 'linear', 'EndValues', 'nearest');
    xB = xB - nanmean(xB);
    dt_bpr_days = days(median(diff(time_list{i}(maskTimeB))));
    if ~isfinite(dt_bpr_days) || dt_bpr_days<=0, dt_bpr_days = 1/24; end
    [PB, fB] = pmtmPH_ori(xB, dt_bpr_days, nw, 0, numel(xB));
    plot(fB, log10(PB)+off, 'Color', abs_colors(i,:), 'LineStyle', '-', 'LineWidth',1, ...
         'DisplayName',[station_names{i} ' BPR']);

    % === ROMS spectrum ===
    maskTimeR = common_time >= t_start & common_time <= t_end;
    xR = fillmissing(roms_abs{i}(maskTimeR), 'linear', 'EndValues', 'nearest');
    xR = xR - nanmean(xR);
    dt_roms_days = days(hours(2));
    [PR, fR] = pmtmPH_ori(xR, dt_roms_days, nw, 0, numel(xR));
    plot(fR, log10(PR)+off, 'Color', abs_colors(i,:), 'LineStyle', '--', 'LineWidth',1, ...
         'DisplayName',[station_names{i} ' ROMS']);
end

ax3 = gca; ax3.FontSize = 16;
xlim([0 1])
xlabel('Frequency (cycles/day)')
ylabel('log_{10} Power (offset)')
title('Spectral Comparison by Station')
legend('NumColumns',2,'Location','southoutside');
grid on

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
