%% Load the ADCP temperature file
meftem = load("H:\Chapter3\LCH\ADCP\Endeavour_MainEndeavourField_AcousticDopplerCurrentProfiler600kHz_Temperature_20220701T000000Z_20250730T043000Z-NaN_clean_avg15minute.mat");
eastem = load("H:\Chapter3\LCH\ADCP\Endeavour_EndeavourEast_ConductivityTemperatureDepth_Temperature_20220701T000000Z_20250730T043000Z-NaN_clean_avg15minute.mat");
mothra = load("H:\Chapter3\LCH\ADCP\Endeavour_Mothra_BottomPressureRecorder_P-SensorTemperature_20220701T000000Z_20250730T043000Z-NaN_clean_avg15minute.mat");
south = load("H:\Chapter3\LCH\ADCP\Endeavour_EndeavourSouth_BottomPressureRecorder_P-SensorTemperature_20220701T000000Z_20250730T043000Z-NaN_clean_avg15minute.mat");
% Convert MATLAB datenums to datetime *once* and use that for plotting
time = datetime(meftem.data.time, 'ConvertFrom', 'datenum');
time_mo = time;
mef_temp = meftem.data.dat;
east_temp = eastem.data.dat;
mothra_temp = mothra.data.dat;
south_temp = south.data.dat;
%%
exclude_start = datetime(2023, 3, 18);
exclude_end   = datetime(2023, 3, 29);
mask = time_mo < exclude_start | time_mo > exclude_end;

time_mo_filtered = time_mo(mask);
mothra_temp_filtered = mothra_temp(mask);

exclude_start = datetime(2024, 1, 25);
exclude_end   = datetime(2024, 1, 26);
mask = time_mo < exclude_start | time_mo > exclude_end;

time_mo_filtered = time_mo(mask);
mothra_temp_filtered = mothra_temp(mask);
%%
% ------------------------
% Step 1: Define sampling info
% ------------------------
dt_minutes = 15;                                % Sampling interval in minutes
Fs = 1 / (dt_minutes * 60);                     % Sampling frequency in Hz
cutoff_days = 14;                               % Desired cutoff in days
cutoff_freq = 1 / (cutoff_days * 24 * 3600);    % Convert to Hz

% ------------------------
% Step 2: Design lowpass Butterworth filter
% ------------------------
Wn = cutoff_freq / (Fs/2);                      % Normalized cutoff frequency
[b, a] = butter(4, Wn, 'low');                  % 4th-order Butterworth filter

% ------------------------
% Step 3: Apply the filter (zero-phase)
% Replace NaNs using linear interpolation
mef_temp    = fillmissing(mef_temp, 'linear');
east_temp   = fillmissing(east_temp, 'linear');
south_temp  = fillmissing(south_temp, 'linear');

% 1. Remove or replace NaNs and Infs
mothra_temp_filtered(~isfinite(mothra_temp_filtered)) = NaN;

% 2. Fill missing values
% Use 'nearest' to deal with edge NaNs, then 'linear' to smooth gaps
mothra_temp_filtered = fillmissing(mothra_temp_filtered, 'nearest');
mothra_temp_filtered = fillmissing(mothra_temp_filtered, 'linear');

%% Apply lowpass filter
% 
% mothra_temp_lp = filtfilt(b, a, mothra_temp_filtered);
% mef_temp_lp    = filtfilt(b, a, mef_temp);
% east_temp_lp   = filtfilt(b, a, east_temp);
% 
% south_temp_lp  = filtfilt(b, a, south_temp);
% % Fill missing values *after* masking

%% no filter result

mothra_temp_lp = mothra_temp_filtered;
mef_temp_lp    =  mef_temp;
east_temp_lp   =east_temp;

south_temp_lp  = south_temp;
%%
% === Step 1: Define sampling parameters ===
dt_minutes = 15;                                % sampling interval
Fs = 1 / (dt_minutes * 60);                     % Hz
cutoff_days = 14;
cutoff_freq = 1 / (cutoff_days * 24 * 3600);    % Hz

% === Step 2: Design lowpass Butterworth filter ===
Wn = cutoff_freq / (Fs/2);                      % normalized cutoff
[b, a] = butter(4, Wn, 'low');                  % 4th order lowpass


% === Step 4: Interpolate all data onto common time (east reference) ===
time_ref = time_mo;  % Assuming all time vectors are the same
t_common = time_ref;

mef_interp    = interp1(time, mef_temp_lp,    t_common, 'linear', 'extrap');
south_interp  = interp1(time, south_temp_lp,  t_common, 'linear', 'extrap');
mothra_interp = interp1(time_mo_filtered, mothra_temp_lp, t_common, 'linear', 'extrap');
east_interp   = interp1(time, east_temp_lp,   t_common, 'linear', 'extrap');


% === Step 5: Compute difference and demean ===
mef_diff    = mef_interp    - east_interp;
south_diff  = south_interp  - east_interp;
mothra_diff = mothra_interp - east_interp;

mef_diff_demean    = mef_diff    - mean(mef_diff, 'omitnan');
south_diff_demean  = south_diff  - mean(south_diff, 'omitnan');
mothra_diff_demean = mothra_diff - mean(mothra_diff, 'omitnan');

% === Step 6: Plot with vertical offsets ===
figure('Position', [100 100 1200 500]); hold on;

offset_mef    = 0;
offset_south  = 0.03;
offset_mothra = 0.06;

plot(t_common, mef_diff_demean    + offset_mef,    'k-', 'LineWidth', 1.2, 'DisplayName', 'MEF - East');
plot(t_common, south_diff_demean  + offset_south,  'm-', 'LineWidth', 1.2, 'DisplayName', 'South - East');
plot(t_common, mothra_diff_demean + offset_mothra, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Mothra - East');

xlabel('Time'); ylabel('Δ Temperature (°C) + Offset');
title('Demeaned Δ Temperature between Stations (Lowpass Filtered, 14-day)');
legend('Location','southoutside','Orientation','horizontal');
grid on;
xlim([datetime(2022,7,1) datetime(2025,7,1)]);
xtickformat('yyyy-MM-dd');

%%
%%% === Filepaths ===
east_fp  = "C:\Users\River\Downloads\Endeavour_stations\east.mat";
south_fp = "C:\Users\River\Downloads\Endeavour_stations\south.mat";
north_fp = "C:\Users\River\Downloads\Endeavour_stations\north.mat";
mef_fp   = "C:\Users\River\Downloads\Endeavour_stations\mef.mat";
mothra_fp= "C:\Users\River\Downloads\Endeavour_stations\mothra.mat";

% === Load and apply Godin filter ===
[time_e, east_godin]   = d(east_fp);
[time_s, south_godin]  = d(south_fp);
[time_n, north_godin]  = d(north_fp);
[time_m, mef_godin]    = d(mef_fp);
[time_mo, mothra_godin]= d(mothra_fp);
%% === Remove Mothra data between 2023-03-19 and 2023-03-23 ===
% === Ensure time_mo is datetime and mask correctly ===
% === Remove Mothra data between 2023-03-18 and 2023-03-23 ===
exclude_start = datetime(2023, 3, 18);
exclude_end   = datetime(2023, 3, 24);
mask = time_mo < exclude_start | time_mo > exclude_end;

time_mo = time_mo(mask);
mothra_godin = mothra_godin(mask);

%% === Apply additional 14-day lowpass Butterworth filter ===
east_lp   = apply_14day_lowpass(time_e, east_godin);
south_lp  = apply_14day_lowpass(time_s, south_godin);
north_lp  = apply_14day_lowpass(time_n, north_godin);
mef_lp    = apply_14day_lowpass(time_m, mef_godin);
mothra_lp = apply_14day_lowpass(time_mo, mothra_godin);
%% === Create timetables with named variable ===
TT_east   = timetable(time_e, east_lp,   'VariableNames', {'east'});
TT_south  = timetable(time_s, south_lp,  'VariableNames', {'south'});
TT_north  = timetable(time_n, north_lp,  'VariableNames', {'north'});
TT_mef    = timetable(time_m, mef_lp,    'VariableNames', {'mef'});
TT_mothra = timetable(time_mo, mothra_lp,'VariableNames', {'mothra'});

% === Synchronize all timetables ===
TT_all = synchronize(TT_east, TT_south, TT_north, TT_mef, TT_mothra, ...
                     'regular', 'linear', 'TimeStep', hours(1));

% === Extract time and variables ===
time_common = TT_all.Properties.RowTimes;
east   = TT_all.east;
south  = TT_all.south;
north  = TT_all.north;
mef    = TT_all.mef;
mothra = TT_all.mothra;

% === Compute and demean differences ===
diff_south  = south  - east; diff_south  = diff_south  - nanmean(diff_south);
diff_north  = north  - east; diff_north  = diff_north  - nanmean(diff_north);
diff_mef    = mef    - east; diff_mef    = diff_mef    - nanmean(diff_mef);
diff_mothra = mothra - east; diff_mothra = diff_mothra - nanmean(diff_mothra);
offsets = 0:1.5:6;
% === Plot ===
figure;
hold on;
plot(time_common, diff_south*100+offsets(1),  'b', 'DisplayName', 'South - East');
plot(time_common, diff_north*100+offsets(2),  'g', 'DisplayName', 'North - East');
plot(time_common, diff_mef*100+offsets(3),    'r', 'DisplayName', 'MEF - East');
plot(time_common, diff_mothra*100+offsets(4), 'm', 'DisplayName', 'Mothra - East');
xlabel('Time');
ylabel('Demeaned Pressure Difference (hPa)');
title('Demeaned Pressure Differences Relative to East');
legend;
grid on;
xlim([datetime(2022,9,1) datetime(2025,7,1)]);
%% plot together
figure('Position',[100 100 1600 600]);           % taller figure
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ------------------------------------------------------------------------
% TOP PANEL – earthquake counts
% ------------------------------------------------------------------------
ax1 = nexttile(t,1);

%ax1 = subplot(211);
hold(ax1,'on');
plot(ax1, time_common, diff_north*100+offsets(1),  'g-','LineWidth', 1.2, 'DisplayName', 'North - East');
plot(ax1, time_common, diff_mef*100+offsets(2),    'k-','LineWidth', 1.2, 'DisplayName', 'MEF - East');
plot(ax1, time_common, diff_south*100+offsets(3),  'm-', 'LineWidth', 1.2,'DisplayName', 'South - East');
plot(ax1, time_common, diff_mothra*100+offsets(4), 'b-', 'LineWidth', 1.2,'DisplayName', 'Mothra - East');
xlabel(ax1,'Time');
ylabel(ax1,'\Delta Pressure Difference (hPa) + Offset');
title(ax1,'Demeaned \Delta Pressure Relative to East (Lowpass Filtered, 14-day)');
legend(ax1,'Location','southoutside','Orientation','horizontal');
grid(ax1,'on');
xlim(ax1,[datetime(2022,9,1) datetime(2025,7,1)]);


ax2 = nexttile(t,2);

%ax2 = subplot(212);
hold(ax2,'on')
offset_mef    = 0;
offset_south  = 0.03;
offset_mothra = 0.06;
plot(ax2, t_common, mef_diff_demean    + offset_mef,    'k-', 'LineWidth', 1.2, 'DisplayName', 'MEF - East');
plot(ax2, t_common, south_diff_demean  + offset_south,  'm-', 'LineWidth', 1.2, 'DisplayName', 'South - East');
plot(ax2, t_common, mothra_diff_demean + offset_mothra, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Mothra - East');
xlabel(ax2,'Time'); ylabel(ax2,'\Delta Temperature (°C) + Offset');
title(ax2,'Demeaned \Delta Temperature Relative to East (Lowpass Filtered, 14-day)');
legend(ax2,'Location','southoutside','Orientation','horizontal');
grid(ax2,'on');
xlim(ax2,[datetime(2022,9,1) datetime(2025,7,1)]);
xtickformat(ax2,'yyyy-MM-dd');
ylim([-0.06,0.1]);
%
maskStart = datetime(2023,9,16);
maskEnd   = datetime(2024,4,18);
maskCol   = [1 0.6 0.6];      % pink RGB
maskAlpha = 0.25;             % 0 = fully transparent, 1 = opaque

axList = [ax1 ax2];   % all axes you want shaded

% for ax = axList
%     axes(ax);                                   % activate axis
%     yl = ylim(ax);     
%     p  = patch(ax, [maskStart maskEnd maskEnd maskStart], ...
%                     [yl(1)    yl(1)  yl(2) yl(2)], ...
%                     maskCol, 'FaceAlpha', maskAlpha, ...
%                     'EdgeColor','none');
%     % keep the patch in the background so lines/bars stay visible
%     uistack(p,'bottom');
% end



linkaxes([ax1, ax2], 'x');

%%
% === Ensure time alignment between pressure and temperature ===
% Interpolate pressure onto temperature timebase (t_common)
diff_mef_interp    = interp1(time_common, diff_mef,    t_common, 'linear', 'extrap');
diff_south_interp  = interp1(time_common, diff_south,  t_common, 'linear', 'extrap');
diff_mothra_interp = interp1(time_common, diff_mothra, t_common, 'linear', 'extrap');

% === Compute correlation coefficients ===
valid_mef    = isfinite(diff_mef_interp)    & isfinite(mef_diff_demean);
valid_south  = isfinite(diff_south_interp)  & isfinite(south_diff_demean);
valid_mothra = isfinite(diff_mothra_interp) & isfinite(mothra_diff_demean);

R_mef = corr(diff_mef_interp(valid_mef), mef_diff_demean(valid_mef));
R_south = corr(diff_south_interp(valid_south), south_diff_demean(valid_south));
R_mothra = corr(diff_mothra_interp(valid_mothra), mothra_diff_demean(valid_mothra));


% === Display correlation results ===
fprintf('Correlation (MEF ΔP vs ΔT):     %.3f\n', R_mef);
fprintf('Correlation (South ΔP vs ΔT):   %.3f\n', R_south);
fprintf('Correlation (Mothra ΔP vs ΔT):  %.3f\n', R_mothra);

%% === Function to apply 14-day lowpass Butterworth filter ===
function data_lp = apply_14day_lowpass(time_vec, data)
    dt_hours = hours(median(diff(time_vec)));
    fs = 1 / dt_hours;
    cutoff = 1 / (14 * 24);  % cycles per hour for 14-day period
    [b, a] = butter(4, cutoff / (fs/2), 'low');
    data_clean = fillmissing(data, 'linear');
    data_lp = filtfilt(b, a, data_clean);
    data_lp = data_lp - nanmean(data_lp);  % zero-center
end
%%
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



