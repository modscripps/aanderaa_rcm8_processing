function rcm = g_aanderaa_reading(sn,df,cf,lt,last_good_data,magnetic_declination,procfigdir)

% RCM = G_AANDEERA_READING Read and calibrate Aanderaa current meter data
%
%   RCM = f_aanderaa_reading(sn,dd,cf)
%
%   INPUT   sn  - serial number
%           df  - data file incl. path
%           cf  - calibration file
%           lt  - gmt time of last sample (six element vector)
%       
%   OUTPUT  rcm - data structure
%
%   Gunnar Voet, APL - UW - Seattle
%   voet@apl.washington.edu
%
%   Created: 02/09/2014

% magnetic_declination = 10.6;

%% Read and calibrate data

% rcm = struct('spd',{},'dir',{},'prs',{},'tmp',{},'con',{},...
%              'u',{},'v',{},'spd_unit',{},'dir_unit',{},'prs_unit',{},...
%              'tmp_unit',{},'con_unit',{},'mooring_id',sn);

inst_id = sn;

fprintf(1,'\n----\nProcessing SN%1d\n----\n',sn);

d = g_calibrate_rcm(sn,df,cf,lt,last_good_data,procfigdir);

% rcm = d;
rcm.sn = sn;


rcm.spd = d.spd;
rcm.dir = d.dir;
% rcm(1).prs = d.prs;
rcm.tmp = d.tmp;
% rcm(1).con = d.con;
% rcm(1).spd_unit = d.spd_unit;
% rcm(1).dir_unit = d.dir_unit;
% rcm(1).prs_unit = d.prs_unit;
% rcm(1).tmp_unit = d.tmp_unit;
% rcm(1).con_unit = d.con_unit;
rcm.time = d.time;
rcm.ref = d.ref;
% end

%% Correct for magnetic declination
% Magnetic deviation is 9 deg in the area of the Faroe Bank Channel

    rcm.dir = rcm.dir+magnetic_declination;
    k = find(rcm.dir>=360);
    rcm.dir(k) = rcm.dir(k)-360;

%% Change speed to m/s

rcm.spd = rcm.spd./100;
rcm.spd_unit = 'm/s';


%% Calculate u (east) and v (north) components

[rcm.u,rcm.v] = g_speeddir2uv(rcm.spd,rcm.dir);

