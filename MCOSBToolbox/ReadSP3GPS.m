% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [sate] = ReadSP3GPS(ipath,DOYth,sate_mark,year)
    list_obs=dir([ipath '\*.sp3']);
    len=length(list_obs);
    % last, current, and next
    if len<3
       error('need at least three SP3 files !!!');
    end
    % DOY
    current_doy  = [num2str(year,'%04d'),num2str(DOYth,'%03d')];
    for i=1:len
        if ~contains(list_obs(i).name,current_doy)
            continue;
        end
        %last
        pr_obs=[ipath '\' list_obs(i-1).name];
        %current
        cu_obs=[ipath '\' list_obs(i).name];
        %next
        nx_obs=[ipath '\' list_obs(i+1).name];  
        [GPS1, ~, ~] = read_precise_eph(pr_obs);
        [GPS2, ~, ~] = read_precise_eph(cu_obs);
        [GPS3, ~, ~] = read_precise_eph(nx_obs);
        gpsx1 = GPS1.X; gpsy1 = GPS1.Y; gpsz1 = GPS1.Z;
        gpsx2 = GPS2.X; gpsy2 = GPS2.Y; gpsz2 = GPS2.Z;
        gpsx3 = GPS3.X; gpsy3 = GPS3.Y; gpsz3 = GPS3.Z;
        numgps=max([size(gpsx1,2),size(gpsx2,2),size(gpsx3,2)]);
        [sate.gpsx,sate.gpsy,sate.gpsz]=interplotation(numgps,gpsx1,gpsy1,gpsz1,gpsx2,gpsy2,gpsz2,gpsx3,gpsy3,gpsz3);
        % Extension
        gpssate=size(sate_mark.gps,2);
        gpsl=size(sate.gpsx,2);
        if gpsl<gpssate
            sate.gpsx(:,gpsl+1:gpssate)=0;
            sate.gpsy(:,gpsl+1:gpssate)=0;
            sate.gpsz(:,gpsl+1:gpssate)=0;
        end
        break;
    end
end
function [interp_x2,interp_y2,interp_z2]=interplotation(satenum,x1,y1,z1,x2,y2,z2,x3,y3,z3)
    interp_x2=zeros(2880,satenum);
    interp_y2=zeros(2880,satenum);
    interp_z2=zeros(2880,satenum);
    x2=[x1(end-3:end,:);x2;x3(1:5,:)];
    y2=[y1(end-3:end,:);y2;y3(1:5,:)];
    z2=[z1(end-3:end,:);z2;z3(1:5,:)];
    m_t=linspace(-40,2880+40,297);
    for i=1:satenum
        for j=1:288     
            tt=m_t(j:j+9);
            x=x2((j:9+j),i)';y=y2((j:9+j),i)';z=z2((j:9+j),i)';
            t0=linspace(m_t(j+4),m_t(j+5)-1,10);
            interp_x2((10*j-9):10*j,i)=interp_lag(tt,x,t0)';
            interp_y2((10*j-9):10*j,i)=interp_lag(tt,y,t0)';
            interp_z2((10*j-9):10*j,i)=interp_lag(tt,z,t0)';
        end
    end
end
function y0 = interp_lag (x, y, x0)
    n=length(x);
    y0=zeros(size(x0));
    for k=1:n
        t=1;
        for i=1:n
            if i~=k
                t=t.*(x0-x(i))/(x(k)-x(i));
            end
        end
        y0=y0+t*y(k);
    end
end
%% ----------------subfunction-------------------
% From the raPPPid by M.F. Glaner
function [GPS, GAL, BDS] = read_precise_eph(filename)
% INPUT:
% 	filename        string with path and name of .sp3 file
% OUPUT:
% 	GPS:          	struct with following elements where
%                       -rows: epochs
%                       -columns: prn numbers of satellites
%       GPS.t:         	time in sec. of week
%       GPS.x:          x coordinate of sat. orbit [m]
%       GPS.y:          y coordinate of sat. orbit [m]
%       GPS.z:          z coordinate of sat. orbit [m]
%       GPS.dt:         clock correction [s]
dbstop if error
GPS = []; 
GAL = []; 
BDS = [];
% open and read sp3-file
fid = fopen(filename);
lines = textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
lines = lines{1};
fclose(fid);
no_lines = length(lines); % number of lines of file

i = 1;
while 1   % loop to go to the first entry
    tline = lines{i};
    i = i + 1;
    if (tline(1) == '*')
        break
    end
end

idx = 0; 	% index of epoch
while i <= no_lines         % loop till end of file   
    if tline(1) ~= '*'
        continue
    end
    epochheader = sscanf(tline,'%*c %f %f %f %f %f %f');      % start with epoch header (always in GPS time)
    date = epochheader;
    jd = cal2jd_GT(date(1), date(2), date(3) + date(4)/24 + date(5)/1440 + date(6)/86400);
    [~, sow,~] = jd2gps_GT(jd);
    idx = idx + 1;          % increase epoch index
       
    tline = lines{i};   
    i = i + 1;
    % loop over data entry
    while tline(1) ~= '*' && tline(1) ~= 'E'            
        % - jump over other GNSS
        if ~contains('GRECJ', tline(2))
            tline = lines{i};   i = i + 1;
            continue
        end
        
        % - get data
        type = tline(1);                % Position or Velocity data
        Epoch = sscanf(tline(3:end),'%f');
        prn = Epoch(1);
        X = Epoch(2);       	
        Y = Epoch(3);
        Z = Epoch(4);
        dT = Epoch(5)*10^-6;        % [microsec] to [s]
        if dT >= 0.1                 % sometimes dt is denoted as 99999.99
            dT = 0;
        end
        gnss = tline(2);
        
        % - save data
        if type == 'P'          % Position data
            if gnss == 'G'
                GPS = save_position(GPS, idx, prn, sow, X, Y, Z, dT);
            elseif gnss == 'E'
                GAL = save_position(GAL, idx, prn, sow, X, Y, Z, dT);
            elseif gnss == 'C'
                BDS = save_position(BDS, idx, prn, sow, X, Y, Z, dT);                
            end
        end    
        % - get next line
        tline = lines{i};   
        i = i + 1;
    end
end
% check read data to prevent errors later-on
% e.g., size of struct and GPS week rollover 
% GPS = checkPrecEph(GPS, 32, bool_check);
% GAL = checkPrecEph(GAL, 62, bool_check);
% BDS = checkPrecEph(BDS, 36, bool_check);

% save-function for position
function GNSS = save_position(GNSS, i, prn, sow, X, Y, Z, dT)
GNSS.t (i, prn)	= round(sow);
GNSS.X (i, prn)	= X*1000;       % [km] to [m]
GNSS.Y (i, prn)	= Y*1000;
GNSS.Z (i, prn) = Z*1000;
GNSS.dT(i, prn) = dT;
end

% check read-in precise ephemeris
function GNSS = checkPrecEph(GNSS, noSats, bool_check)
if isempty(GNSS)
    return              % no data for this GNSS
end
[rows, sats] = size(GNSS.t);
% check for missing columns/satellites
if sats < noSats       
    GNSS.t (rows, noSats) = 0;
    GNSS.X (rows, noSats) = 0;
    GNSS.Y (rows, noSats) = 0;
    GNSS.Z (rows, noSats) = 0;
    GNSS.dT(rows, noSats) = 0;
end
% check for GPS week rollover:
% If the processed day is a Saturday, the last epoch of the precise orbit 
% file is 0h on the next day = Sunday = 1st day of new GPS week. Therefore 
% the time-stamp of this epoch has to be corrected from 0 to 604800 [sow]
newGPSweek = GNSS.t < GNSS.t(1,:);
GNSS.t = GNSS.t + newGPSweek.*86400*7;
if bool_check
% check if there are epochs without data at all, these epochs prevent 
% reasonable interpolation during processing -> delete these epochs
    bool_empty_epoch =  all(GNSS.X == 0, 2)  | all(GNSS.Y == 0, 2)  | all(GNSS.Z == 0, 2) | all(isnan(GNSS.X),2) | all(isnan(GNSS.Y),2) | all(isnan(GNSS.Z),2);
    GNSS.t(bool_empty_epoch,:) = [];
    GNSS.X(bool_empty_epoch,:) = [];
    GNSS.Y(bool_empty_epoch,:) = [];
    GNSS.Z(bool_empty_epoch,:) = [];
    GNSS.dT(bool_empty_epoch,:) = [];
end
end
end