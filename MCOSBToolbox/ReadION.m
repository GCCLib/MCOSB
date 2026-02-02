% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [ION] = ReadION(ipath,DOYth)
    list_1=dir([ipath '/*.*i']);
    list_2=dir([ipath '/*.*inx']);
    list=[list_1;list_2];
    len=length(list);
    if len<1
       error('need at least a ion files !!!');
    end
    % DOY
    current_doy  = [num2str(DOYth,'%03d'),'0'];
    for i=1:len
        if ~contains(list(i).name,current_doy)
            continue;
        end
        path_ionex = [ipath '/' list(i).name];
        ION.ionex = read_ionex_TUW(path_ionex);
        break;
    end
end

function [IonexFile] = read_ionex_TUW(filename)
% Reading in IONEX files and save TEC maps and header information to struct variable. RMS-Maps are not read in. The time system of IONEX maps is UTC
% time, details on the format: https://files.igs.org/pub/data/format/ionex1.pdf based on function created on 15.11.2017 by Janina Boisits
% INPUT:
%   filename        string, path to IONEX file
% OUTPUT:
%   IonexFile       struct, contains data of IONEX file
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner

% read file
fid   = fopen(filename);
ionex = textscan(fid,'%s','delimiter', '\n', 'whitespace', '');
fclose(fid);
ionex = ionex{1};
ionex_char = char(ionex);

% extract start and end epoch
idx = contains(ionex,'EPOCH OF FIRST MAP');
start_epoch(1) = sscanf(ionex{idx}(1:6), '%f');
start_epoch(2) = sscanf(ionex{idx}(7:12), '%f');
start_epoch(3) = sscanf(ionex{idx}(13:18), '%f');
start_epoch(4) = sscanf(ionex{idx}(19:24), '%f');
start_epoch(5) = sscanf(ionex{idx}(25:30), '%f');
start_epoch(6) = sscanf(ionex{idx}(31:36), '%f');
[~,start_sow,~] = jd2gps_GT(cal2jd_GT(start_epoch(1), start_epoch(2), start_epoch(3) + start_epoch(4)/24 + start_epoch(5)/1440 + start_epoch(6)/86400));
idx = contains(ionex,'EPOCH OF LAST MAP');
end_epoch(1)   = sscanf(ionex{idx}(1:6), '%f');
end_epoch(2)   = sscanf(ionex{idx}(7:12), '%f');
end_epoch(3)   = sscanf(ionex{idx}(13:18), '%f');
end_epoch(4)   = sscanf(ionex{idx}(19:24), '%f');
end_epoch(5)   = sscanf(ionex{idx}(25:30), '%f');
end_epoch(6)   = sscanf(ionex{idx}(31:36), '%f');
[~,end_sow,~]  = jd2gps_GT(cal2jd_GT(end_epoch(1), end_epoch(2), end_epoch(3) + end_epoch(4)/24 + end_epoch(5)/1440 + end_epoch(6)/86400));
if end_sow == 0 && end_sow < start_sow
    % set saved 'end of ionex file' to the last second of the GPS week to 
	% avoid errors later on (e.g., in iono_gims.m)
    end_sow = 604800;       
end
% extract other header information
idx = contains(ionex,'INTERVAL');
Interval      = sscanf(ionex{idx}(1:6), '%f');
idx = contains(ionex,'# OF MAPS IN FILE');
mapsNR        = sscanf(ionex{idx}(1:6), '%f');
idx = contains(ionex,'MAPPING FUNCTION');
MappingFun    = ionex{idx}(1:60);
MappingFun(MappingFun == ' ') = '';         % remove empty spaces from mapping function
idx = contains(ionex,'ELEVATION CUTOFF');
CutOff_EL     = sscanf(ionex{idx}(1:8), '%f');
idx = contains(ionex,'BASE RADIUS');
EarthRadius   = sscanf(ionex{idx}(1:8), '%f');
idx = contains(ionex,'MAP DIMENSION');
MapDimension  = sscanf(ionex{idx}(1:6), '%f');
idx = contains(ionex,'HGT1 / HGT2 / DHGT');
hgt1          = sscanf(ionex{idx}(1:8), '%f');
hgt2          = sscanf(ionex{idx}(9:14), '%f');
dhgt          = sscanf(ionex{idx}(15:20), '%f');
idx = contains(ionex,'LAT1 / LAT2 / DLAT');
lat1          = sscanf(ionex{idx}(1:8), '%f');
lat2          = sscanf(ionex{idx}(9:14), '%f');
dlat          = sscanf(ionex{idx}(15:20), '%f');
idx = contains(ionex,'LON1 / LON2 / DLON');
lon1          = sscanf(ionex{idx}(1:8), '%f');
lon2          = sscanf(ionex{idx}(9:14), '%f');
dlon          = sscanf(ionex{idx}(15:20), '%f');
idx = contains(ionex,'EXPONENT');
exponent      = sscanf(ionex{idx}(1:6), '%f');

% dimension of TEC maps
latNR   = length(lat1:dlat:lat2);
lonNR   = length(lon1:dlon:lon2);
maps    = NaN(latNR,lonNR,mapsNR);
linesNR = ceil(lonNR/16);

% read TEC maps
idxStart = find(contains(ionex,'START OF TEC MAP'));
idxEnd   = find(contains(ionex,'END OF TEC MAP'));
for j = 1:mapsNR
    currMap      = ionex_char(idxStart(j):idxEnd(j),:);
    currMap_cell = cellstr(currMap);
    idx          = find(contains(currMap_cell,'LAT/LON1/LON2/DLON/H'));
    for ilat = 1:latNR
        latBand        = currMap(idx(ilat)+1:idx(ilat)+linesNR,:);
        latValues      = textscan(latBand','%f');
        if numel(maps(ilat,:,j)) ~= numel(latValues{1}) 	% TEC values > 10^4
            [rows, cols] = size(latBand);
            empty = latBand(:,1); empty(:) = ' ';           % creaty empty block
            pos = 1;    % loop to insert whitespaces to seperate TEC values > 10^4
            while pos < size(latBand,2)     
                latBand = [latBand(:,1:pos+4), empty, latBand(:,pos+5:end)];
                pos = pos + 6;
            end
            latValues      = textscan(latBand','%f');       % reread with seperated values
        end
        maps(ilat,:,j) = latValues{1};
    end
end
% safe ionex file to struct
IonexFile = struct('start_epoch',start_epoch, 'end_epoch',end_epoch,...
    'start_sow',start_sow, 'end_sow',end_sow,...
    'interval',Interval, 'number_of_maps',mapsNR, 'mf',MappingFun,...
    'cutoff',CutOff_EL, 'radius',EarthRadius,...
    'map_dimension',MapDimension, 'hgt',[hgt1,hgt2,dhgt],...
    'lat',[lat1,lat2,dlat], 'lon',[lon1,lon2,dlon],...
    'exponent',exponent, 'map',maps);
end

function jd = cal2jd_GT(yr, mn, dy)
    % cal2jd_GT  Converts calendar date to Julian date using algorithm from "Practical Ephemeris Calculations" by Oliver Montenbruck
    %   (Springer-Verlag, 1989). Uses astronomical year for B.C. dates
    %   (2 BC = -1 yr). Non-vectorized version. See also doy2jd_GT, gps2jd_GT,
    %   jd2cal_GT, jd2dow_GT, jd2doy_GT, jd2gps_GT, jd2yr_GT, yr2jd_GT.
    % Usage:   jd=cal2jd_GT(yr,mn,dy)
    % Input:   yr - calendar year (4-digit including century)
    %          mn - calendar month
    %          dy - calendar day (including fractional day)
    % Output:  jd - jJulian date  
    % Copyright (c) 2011, Michael R. Craymer
    % All rights reserved.
    % Email: mike@craymer.com
    
    jd = 0;   
    if nargin ~= 3
      warning('Incorrect number of input arguments');
      return;
    end
    if mn < 1 | mn > 12
      warning('Invalid input month');
      return
    end
    if dy < 1
      if (mn == 2 & dy > 29) | (any(mn == [3 5 9 11]) & dy > 30) | (dy > 31)
        warning('Invalid input day');
        return
      end
    end   
    if mn > 2
      y = yr;
      m = mn;
    else
      y = yr - 1;
      m = mn + 12;
    end
    
    date1=4.5+31*(10+12*1582);   % Last day of Julian calendar (1582.10.04 Noon)
    date2=15.5+31*(10+12*1582);  % First day of Gregorian calendar (1582.10.15 Noon)
    date=dy+31*(mn+12*yr);    
    if date <= date1
      b = -2;
    elseif date >= date2
      b = fix(y/400) - fix(y/100);
    else
      warning('Dates between October 5 & 15, 1582 do not exist');
      return;
    end    
    if y > 0
      jd = fix(365.25*y) + fix(30.6001*(m+1)) + b + 1720996.5 + dy;
    else
      jd = fix(365.25*y-0.75) + fix(30.6001*(m+1)) + b + 1720996.5 + dy;
    end
end