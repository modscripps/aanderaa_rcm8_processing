% TTIDE
%
% Post-cruise processing of Aanderaa RCM data
%
% Gunnar Voet
% gvoet@ucsd.edu
%
% Created: 09/07/2015

clear

% Directory for cal-files
cfdir = 'cal_files/';

addpath functions/
addpath more_functions/

sn = 4917;
mn = 'M5';
cf = sprintf('%scal_%1d',cfdir,sn);
% df = '/Users/gunnar/Projects/ttide/data/Moorings/M5/Aanderaa/DSU4917_M5_TTIDE2015_corrected.Asc';
% df = '/Volumes/Ahua/data_archive/WaveChasers-DataArchive/TTIDE/Moorings/M5/Aanderaa/DSU4917_M5_TTIDE2015_Corrected.Asc';
df = 'raw/DSU4917_M5_TTIDE2015_Corrected.Asc';
lt = [2015 3 7 2 39 0];

%% Read Aanderaa ascii file

% Write screen output to log file
LogName = sprintf('log/SN%1d_%s.txt',sn,mn);  % log file name
if exist(LogName,'file')==2
system(sprintf('rm %s',LogName));             % delete prev log file
end
system(sprintf('touch %s',LogName));          % create empty log file
diary(LogName)                                % turn logging on

% Read and calibrate data
magnetic_declination = (15+12./60);
rcm = g_aanderaa_reading(sn,df,cf,lt,NaN,magnetic_declination,'procfig');

diary off                                     % turn logging off

% Add mooring name to structure
rcm.mn = mn;

% cut structure
rcm = g_cutstruct(rcm,[677:15354]);

%% Save to .mat file
SaveName = sprintf('proc/SN%1d_%s.mat',sn,mn);
save(SaveName,'rcm')


%% Plot

figure(1)
clf
subaxis(4,1,1)
plot(rcm.time,rcm.u)
grid on
ylabel('u [m/s]')
tlabel
subaxis(4,1,2)
plot(rcm.time,rcm.v)
grid on
ylabel('v [m/s]')
tlabel
subaxis(4,1,3)
plot(rcm.time,rcm.spd)
grid on
ylabel('speed [m/s]')
tlabel
subaxis(4,1,4)
plot(rcm.time,rcm.dir)
grid on
ylabel('direction [deg]')
tlabel