function rcm = g_calibrate_rcm(sn,df,cf,lt,last_good_data,procfigdir)

% rcm = g_calibrate_rcm(inst_id)
%
% Function for the readout and calibration of Aanderaa RCM data. Needs the
% function g_read_rcm and a file cal_<inst_id> with the calibration
% coefficients.
%
% Input: instrument id (numeric).
%        filename - name of the ascii input file.
%
% Output: rcm - structured array with calibrated RCM data.
%
% by G.Voet, IfM Hamburg
% gunnar.voet@zmaw.de
%
% last modification: 05.08.2008

PlotFigure = 1;

%% Read the calibration coefficients

run(cf);

%% Read the ascii data
d = g_read_rcm(df);

length_of_time_series = length(d(:,1));
fprintf(1,'\nTime series has %1d data points.\n\n',...
          length_of_time_series)

% Cut time series if last good data point is given
d1 = d;
if isfinite(last_good_data)
  d = d(1:last_good_data,:);
  fprintf(1,'Cutting time series from data point %1d to %1d\n\n',...
            last_good_data,length_of_time_series);
end


%% Calibrate data and assign to structured array
if rcmcal.id==sn;
% Reference

c = rcmcal.ref_channel;
rcm.ref =   rcmcal.ref(1) + ...
            rcmcal.ref(2).*d(:,c) + ...
            rcmcal.ref(3).*d(:,c).^2 + ...
            rcmcal.ref(4).*d(:,c).^3;
rcm.ref =   rcm.ref';

kref = find(rcm.ref < rcmcal.ref_reading-1 | ...
            rcm.ref > rcmcal.ref_reading+1);

figure(1)
clf
subaxis(4,1,1)
plot(rcm.ref,'k')
hold on
text(0.1,0.9,sprintf('ref reading: %1d',rcmcal.ref_reading),...
     'units','normalized')
set(gca,'ylim',[rcmcal.ref_reading-5 rcmcal.ref_reading+5])


% Temperature
if rcmcal.tmp_channel>0
    c = rcmcal.tmp_channel;
    rcm.tmp =   rcmcal.tmp(1) + ...
                rcmcal.tmp(2).*d(:,c) + ...
                rcmcal.tmp(3).*d(:,c).^2 + ...
                rcmcal.tmp(4).*d(:,c).^3;
    rcm.tmp =   rcm.tmp';
    rcm.tmp_unit = rcmcal.tmp_unit;
    
    rcm.tmp(kref) = NaN;
    
    k = find(d(:,c)==1023);
    rcm.tmp(k) = NaN;
    
else
    rcm.tmp = [];
    rcm.tmp_unit = [];
end

% Pressure
if rcmcal.prs_channel>0
    c = rcmcal.prs_channel;
    rcm.prs =   rcmcal.prs(1) + ...
                rcmcal.prs(2).*d(:,c) + ...
                rcmcal.prs(3).*d(:,c).^2 + ...
                rcmcal.prs(4).*d(:,c).^3;
    rcm.prs =   rcm.prs';
    rcm.prs_unit = rcmcal.prs_unit;
    
    rcm.prs(kref) = NaN;
    
    % Convert MPa into dbar
    if strcmp(rcm.prs_unit,'MPa')
        rcm.prs = rcm.prs.*100;
        rcm.prs_unit = 'dbar';
    elseif strcmp(rcm.prs_unit,'kg/cm2')
        rcm.prs = rcm.prs.*10;
        rcm.prs_unit = 'dbar';
    end
    
else
    rcm.prs = [];
    rcm.prs_unit = [];
end

% Conductivity
if rcmcal.con_channel>0
    c = rcmcal.con_channel;
    rcm.con =   rcmcal.con(1) + ...
                rcmcal.con(2).*d(:,c) + ...
                rcmcal.con(3).*d(:,c).^2 + ...
                rcmcal.con(4).*d(:,c).^3;
    rcm.con =   rcm.con';
    rcm.con_unit = rcmcal.con_unit;
    
    rcm.con(kref) = NaN;
    
    % Convert mmho/cm to S/m
    if strcmp(rcm.con_unit,'mmho/cm')
        rcm.con = rcm.con./10;
        rcm.con_unit = 'S/m';
    end
    
else
    rcm.con = [];
    rcm.con_unit = [];
end

% Speed
if rcmcal.spd_channel>0
    c = rcmcal.spd_channel;
    rcm.spd =   rcmcal.spd(1) + ...
                rcmcal.spd(2).*d(:,c) + ...
                rcmcal.spd(3).*d(:,c).^2 + ...
                rcmcal.spd(4).*d(:,c).^3;
    rcm.spd =   rcm.spd';
    rcm.spd_unit = rcmcal.spd_unit;
    
    rcm.spd(kref) = NaN;
    
else
    rcm.spd = [];
    rcm.spd_unit = [];
end

% Direction
if rcmcal.dir_channel>0
    c = rcmcal.dir_channel;
    rcm.dir =   rcmcal.dir(1) + ...
                rcmcal.dir(2).*d(:,c) + ...
                rcmcal.dir(3).*d(:,c).^2 + ...
                rcmcal.dir(4).*d(:,c).^3;
    rcm.dir =   rcm.dir';
    rcm.dir_unit = rcmcal.dir_unit;
    
    rcm.dir(kref) = NaN;
    
else
    rcm.dir = [];
    rcm.dir_unit = [];
end


%% Inspect low speeds


Sraw = rcm.spd;


% CURRENTS ----------
% Handle near-zero current speeds
% Examining S reveals discrete lower values
Slowervalues = zeros(10,1);
SS = Sraw;
for i = 1:10
	Slowervalues(i) = min(SS);
	j = find(SS(:)==Slowervalues(i));
	SS(j) = NaN;
end
Slowervalues
Sdifferences = Slowervalues(2:10) - Slowervalues(1:9)

S = Sraw;
Smin = min(Sraw);
Sminhalf = Smin/2;
%j = find(Sraw(:)==Smin);
j = find(Sraw(:)<=Smin);
n = length(j);

k = [];
for i = 2:n-1
	if (j(i-1)+1 == j(i)) & (j(i)+1 == j(i+1))
		k = [k; j(i)];
	end
end
S(k) = Sminhalf;
% save SminStudy_NP08_6543 SminStudy_NP08_6543 -ascii

% Get spans of near-zero current
disp('Spans of near-zero current')
kspans = [];
i = 1;
nk = length(k)
k1 = k(1:nk-1);
k2 = k(2:nk);
j = find(k1(:)+1 ~= k2(:));
kspans = [[1;k2(j)], [k1(j);k(nk)]];
kspans = [kspans, kspans(:,2)-kspans(:,1)+1]

subaxis(4,1,2)
plot(S,'r')
hold on
plot(Sraw,'k');
set(gca,'ylim',[0 20])

rcm.spd = S;


%% The time vector
time_temp = d(:,end)';

k = ~isnan(time_temp);
k2 = find(k);

subaxis(4,1,3)
plot(diff(time_temp(k)).*24-nanmedian(diff(time_temp(k))).*24);
ylabel('\Delta t anomaly [h]')
set(gca,'ylim',[-3 3])

rcm.sampling_freq = 288;

% % Do not correct if time difference is less than 1 hour, simply expand the
% % time vector.
% tlength = length(rcm.spd)-1/(rcm.sampling_freq);
% rcm.time = time_temp(1):1/(rcm.sampling_freq):time_temp(1)+tlength/(rcm.sampling_freq);
% % Difference between expanded time vector and last time stamp
% deltat = time_temp(k2(end))-rcm.time(k2(end));
% if abs(deltat) > 1/24
%     error('Time difference too big!')
% end

% Better: Interpolate between the time stamps. There are a few jumps in time
% that correspond with missing (or too much) data, so interpolation is
% better!
rcm.time2 = nan(size(time_temp));
for i = 1:length(k2)-1
  ka = k2(i);
  kb = k2(i+1);
  x = ka:kb;
  rcm.time2(x) = interp1([ka kb],time_temp([ka kb]),x);
end
% Fill last day
% nld = length(rcm.time)-k2(end);
xend = k2(end)+1:length(time_temp);
for i = 1:length(xend)
rcm.time2(xend(i)) = time_temp(k2(end))+i*(1/rcm.sampling_freq);
end


% ok, use the interpolated time
rcm.time = rcm.time2;
rcm = rmfield(rcm,'time2');


fprintf(1,'Readout and calibration of instrument %1d successful!\n\n',sn);

else
    display('Wrong calibration data!')
end

