function rcm = g_read_rcm(dateiname)

% rcm = g_read_rcm(dateiname)
%
% Funktion zum Auslesen von RCM Ascii-Daten.
%
% Input: dateiname
%
% Output: rcm: Struktur mit unkalibrierten RCM Daten.
%
% Teile uebernommen von J.Holfort - RCMlese.m.
%
% by G.Voet, IfM Hamburg
% gunnar.voet@zmaw.de
%
% last modification: 04.08.2008


rcma = [];
fid = fopen(dateiname);
l=' ';
idata=0;
jdata=0;
    
if fid<=0; l=1; end
while ischar(l);
      l = fgetl(fid);
      jfehler=0;
      eval('d1=str2num(l);','jfehler=1;'); % Das scheint so konstruiert zu
                                           % sein, dass der zweite Befehl
                                           % ausgefuehrt wird wenn der
                                           % erste nicht funktioniert.
    if jfehler==0
        idata = idata+1;
        if jdata < length(d1)              % Die erste Zeile ist eine Zeile
                                           % mit Zeitinformation.
            for ii = 1:length(d1)
                rcma(idata,ii) = d1(ii);
            end
            jdata = length(d1);
        else
            rcma(idata,1:length(d1)) = d1;
        end
    end
end
if fid>0;
    fclose(fid);
end

%% search for time stamps and append time, then delete time stamps
jtime=0;

kk = ~(rcma(:,1) == 7);
if jdata >= 6 && sum(~kk) > 1
    rcma = [rcma zeros(idata,1)+NaN];
    jtime = jdata+1;
    for i = 1:idata-1
        if rcma(i,1) == 7;
            jahr = rcma(i,2)+1900;
            if jahr < 1970;
                jahr = jahr+100;
            end
            datum = datenum([jahr,rcma(i,3),rcma(i,4),rcma(i,5),rcma(i,6),0]);
            rcma(i+1,jdata+1) = datum;
        end
    end
else
display('no time decoding posible')
end
rcm = rcma(kk,:);
% idata = sum(kk);
