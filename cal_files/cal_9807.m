% Calibration coefficients for Aanderaa RCM.
%
% Mooring M4
%
% SN 9807
%
% Set channels to zero if sensor not existent.
%
% by G.Voet, IfM Hamburg
% gunnar.voet@zmaw.de
%
% last modification: 09/02/2014

rcmcal.id = 9807;

rcmcal.ref_reading = 783;
rcmcal.ref = [0 1 0 0];
rcmcal.ref_channel = 1;

rcmcal.tmp = [-2.749 8.195e-3 -1.601e-7 7.991e-11];
rcmcal.tmp_channel = 4;
rcmcal.tmp_unit = 'deg C';

rcmcal.con = [0 0 0 0];
rcmcal.con_channel = 0;
rcmcal.con_unit = 'mmho/cm';

rcmcal.prs = [0 0 0 0];
rcmcal.prs_channel = 0;
rcmcal.prs_unit = 'MPa';

rcmcal.dir = [1 3.5e-1 0 0];
rcmcal.dir_channel = 5;
rcmcal.dir_unit = 'Deg. M';

rcmcal.spd = [1.1 2.906e-1 0 0];
rcmcal.spd_channel = 6;
rcmcal.spd_unit = 'cm/s';