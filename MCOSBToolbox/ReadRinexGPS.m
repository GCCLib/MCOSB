% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function ReadRinexGPS(ipath,opath,STA_opath)
    list_obs1=dir([ipath '/*.rnx']);
    list_obs2=dir([ipath '/*.*o']);
    list_obs=[list_obs1;list_obs2];
    if isempty(list_obs)
        sprintf('No RINEX files found in %s',ipath);
        return;
    end
    len=length(list_obs);
    if exist(opath,'dir')==0 
        mkdir(opath);
    end
    if exist(STA_opath,'dir')==0 
        mkdir(STA_opath);
    end
    % Station information
    STAInformation.name=cell(1,len);
    STAInformation.doy=linspace(0,0,len)';
    STAInformation.coor=zeros(len,3);
    STAInformation.RCVType=cell(1,len);
    %--------------------------------------------------------------------------
    for i=1:len
        obsn=list_obs(i).name;
        % 15-s obs file
        if contains(obsn,'15S')
            continue;
        end
        % empty
        if list_obs(i).bytes == 0
            continue;
        end
        % Long file
        if length(obsn)>25
            doy=obsn(17:19);
            STAInformation.name{1,i}=upper(obsn(1:4));
            % OBS
            name1=[upper(obsn(1:4)) doy '.mat'];
            if exist([opath,'/',doy],'dir')==0
                mkdir([opath,'/',doy]);
            end
            filename1=[opath,'/',doy,'/',name1];
            % STA
            name2=['STAInformation' doy '.mat'];
            if exist([STA_opath,'/',doy],'dir')==0
                mkdir([STA_opath,'/',doy]);
            end
            filename2=[STA_opath,'/',doy,'/',name2];
        % Short file
        else
            doy=[obsn(10:11),obsn(5:7)];
            STAInformation.name{1,i}=upper(obsn(1:4));
            % OBS
            name1=[upper(obsn(1:4)) doy '.mat'];
            if exist([opath,'/',doy],'dir')==0
                mkdir([opath,'/',doy]);
            end
            filename1=[opath,'/',doy,'/',name1];
            % STA
            name2=['STAInformation' doy '.mat'];
            if exist([STA_opath,'/',doy],'dir')==0
                mkdir([STA_opath,'/',doy]);
            end
            filename2=[STA_opath,'/',doy,'/',name2];
        end
        % DOY
        STAInformation.doy(i)=str2double(doy);
        % ---------------------------------------------------------------------
        fname=[ipath '/' obsn];
        [obs,coor,RCVType]=read_rnx3GPS(fname); 
        % Save OBS
        if  ~isempty(obs) && sum(coor) ~= 0
            save(filename1,'obs','coor','-mat');
            STAInformation.coor(i,:)=coor;
            STAInformation.RCVType{1,i}=RCVType;
        end
        disp(['-------- [',num2str(doy),'] [ ',num2str(i),'/',num2str(len),' ] GPSOBS file is read !']);
    end
    save(filename2,'STAInformation','-mat');
end
% Read OBS from RINEX3 files
function [obs,coor,RCVType]=read_rnx3GPS(path)
% INPUT: 
%   path: path of observation files
% OUTPUT: 
%   obs:  struct of observation files
%   coor: station coordinates
%   RCVType: receiver type
GPSFlag = 0;
gpsmark = 0;
fid=fopen(path,'r');
while 1
    line=fgetl(fid);
    if ~ischar(line)
        break;
    end  
    % the rinex file is null or end
    % RINEX version should be larger than 3.02
    if length(line)>79 && strcmp(line(61:80),'RINEX VERSION / TYPE')
        if strcmp(line(6),'3')==0
            obs=[];coor=[0,0,0]; RCVType='';
            break;
        end
        if str2double(line(9))<3
            obs=[];coor=[0,0,0];RCVType='';
            break;
        end
    end
    % get receivers coordinates
    if length(line)>78 && strcmp(line(61:79),'APPROX POSITION XYZ')
        coor(1)=str2double(line(1:14));  % X (WGS-84)
        coor(2)=str2double(line(15:28)); % Y (WGS-84)
        coor(3)=str2double(line(29:42)); % Z (WGS-84)
        continue;
    end
    % get receiver types
    if length(line)>78 && strcmp(line(61:79),'REC # / TYPE / VERS')
        RCVType = line(21:32);
        continue;
    end
    %  get GPS observable types
    if length(line)>78 && strcmp('G',line(1)) && strcmp('SYS / # / OBS TYPES',line(61:79))
        gpsloc=zeros(1,18);
        obsgps_n=str2double(line(4:6));
        if obsgps_n>39
            for i=1:13
                gpst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                gpst(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                gpst(i+26)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsgps_n-39
                gpst(i+39)={line(4+4*i:6+4*i)};
            end
        elseif obsgps_n>26 && obsgps_n<=39
            for i=1:13
                gpst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                gpst(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsgps_n-26
                gpst(i+26)={line(4+4*i:6+4*i)};
            end
        elseif obsgps_n>13 && obsgps_n<=26
            for i=1:13
                gpst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsgps_n-13
                gpst(i+13)={line(4+4*i:6+4*i)};
            end
        else
            for i=1:obsgps_n
                gpst(i)={line(4+4*i:6+4*i)};
            end
        end
        % L1：C1C/C1S/C1L/C1X/C1P/C1W/C1N
        % L2：C2C/C2D/C2S/C2L/C2X/C2P/C2W/C2N
        % L5：C5I/C5Q/C5X
        for i=1:obsgps_n
            if strcmp(gpst(1,i),'C1C')
                gpsloc(1)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C1S')
                gpsloc(2)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C1L')
                gpsloc(3)=i;
                continue;
            end            
            if strcmp(gpst(1,i),'C1X')
                gpsloc(4)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C1P') 
                gpsloc(5)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C1W')
                gpsloc(6)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C1N')
                gpsloc(7)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C2C')
                gpsloc(8)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C2D')
                gpsloc(9)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C2S')
                gpsloc(10)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C2L')
                gpsloc(11)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C2X')
                gpsloc(12)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C2P')
                gpsloc(13)=i;
                continue;
            end            
            if strcmp(gpst(1,i),'C2W')
                gpsloc(14)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C2N')
                gpsloc(15)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C5I')
                gpsloc(16)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C5Q')
                gpsloc(17)=i;
                continue;
            end
            if strcmp(gpst(1,i),'C5X')
                gpsloc(18)=i;
                continue;
            end
        end
        % Note: 
        % the following code types should be observed for GPS
        % at least C1W/C2W/C5Q, C1W/C2W/C5X, or C1W/C2W/C5I
        if (gpsloc(6)>0&&gpsloc(14)>0&&gpsloc(17)>0) || (gpsloc(6)>0&&gpsloc(14)>0&&gpsloc(18)>0) || (gpsloc(6)>0&&gpsloc(14)>0&&gpsloc(16)>0)
            GPSFlag = 1;
        end
        continue;
    end      
    % start getting observables
    if(length(line)>72&&strcmp(line(61:73),'END OF HEADER')) 
        % Triple-frequency observations are not existed
        if GPSFlag == 0
            obs=[]; coor=[0,0,0]; RCVType='';
            break;
        end
        
        while 1                      
            line=fgetl(fid);
            if ~ischar(line)
                % come to end,break
                break;
            end 
            % get epoch number
            if length(line)>32 && strcmp(line(1),'>')
                h=str2double(line(14:15));
                m=str2double(line(17:18));
                s=str2double(line(20:29));
                % 30-s sampling rate
                ep=h*120+m*2+s/30+1; 
            end 
            if ep~=fix(ep),continue,end
            if exist('gpsloc','var')
                if strcmp(line(1),'G')
                    gpsmark = 1;
                    sv_G=str2double(line(2:3));
                    if isnan(sv_G) 
                        continue; 
                    end      
                    if gpsloc(1)~=0&&length(line)>gpsloc(1)*16    obs.GPSC1C(ep,sv_G)=str2double(line(16*gpsloc(1)-12:16*gpsloc(1)+1)); end
                    if gpsloc(2)~=0&&length(line)>gpsloc(2)*16    obs.GPSC1S(ep,sv_G)=str2double(line(16*gpsloc(2)-12:16*gpsloc(2)+1)); end
                    if gpsloc(3)~=0&&length(line)>gpsloc(3)*16    obs.GPSC1L(ep,sv_G)=str2double(line(16*gpsloc(3)-12:16*gpsloc(3)+1)); end
                    if gpsloc(4)~=0&&length(line)>gpsloc(4)*16    obs.GPSC1X(ep,sv_G)=str2double(line(16*gpsloc(4)-12:16*gpsloc(4)+1)); end
                    if gpsloc(5)~=0&&length(line)>gpsloc(5)*16    obs.GPSC1P(ep,sv_G)=str2double(line(16*gpsloc(5)-12:16*gpsloc(5)+1)); end
                    if gpsloc(6)~=0&&length(line)>gpsloc(6)*16    obs.GPSC1W(ep,sv_G)=str2double(line(16*gpsloc(6)-12:16*gpsloc(6)+1)); end
                    if gpsloc(7)~=0&&length(line)>gpsloc(7)*16    obs.GPSC1N(ep,sv_G)=str2double(line(16*gpsloc(7)-12:16*gpsloc(7)+1)); end
                    if gpsloc(8)~=0&&length(line)>gpsloc(8)*16    obs.GPSC2C(ep,sv_G)=str2double(line(16*gpsloc(8)-12:16*gpsloc(8)+1)); end
                    if gpsloc(9)~=0&&length(line)>gpsloc(9)*16    obs.GPSC2D(ep,sv_G)=str2double(line(16*gpsloc(9)-12:16*gpsloc(9)+1)); end
                    if gpsloc(10)~=0&&length(line)>gpsloc(10)*16  obs.GPSC2S(ep,sv_G)=str2double(line(16*gpsloc(10)-12:16*gpsloc(10)+1)); end
                    if gpsloc(11)~=0&&length(line)>gpsloc(11)*16  obs.GPSC2L(ep,sv_G)=str2double(line(16*gpsloc(11)-12:16*gpsloc(11)+1)); end
                    if gpsloc(12)~=0&&length(line)>gpsloc(12)*16  obs.GPSC2X(ep,sv_G)=str2double(line(16*gpsloc(12)-12:16*gpsloc(12)+1)); end                   
                    if gpsloc(13)~=0&&length(line)>gpsloc(13)*16  obs.GPSC2P(ep,sv_G)=str2double(line(16*gpsloc(13)-12:16*gpsloc(13)+1)); end
                    if gpsloc(14)~=0&&length(line)>gpsloc(14)*16  obs.GPSC2W(ep,sv_G)=str2double(line(16*gpsloc(14)-12:16*gpsloc(14)+1)); end
                    if gpsloc(15)~=0&&length(line)>gpsloc(15)*16  obs.GPSC2N(ep,sv_G)=str2double(line(16*gpsloc(15)-12:16*gpsloc(15)+1)); end
                    if gpsloc(16)~=0&&length(line)>gpsloc(16)*16  obs.GPSC5I(ep,sv_G)=str2double(line(16*gpsloc(16)-12:16*gpsloc(16)+1)); end
                    if gpsloc(17)~=0&&length(line)>gpsloc(17)*16  obs.GPSC5Q(ep,sv_G)=str2double(line(16*gpsloc(17)-12:16*gpsloc(17)+1)); end
                    if gpsloc(18)~=0&&length(line)>gpsloc(18)*16  obs.GPSC5X(ep,sv_G)=str2double(line(16*gpsloc(18)-12:16*gpsloc(18)+1)); end
                end
            end
        end
        if gpsmark == 0
            obs=[];coor=[0,0,0]; RCVType='';        
        end
    end
end
fclose(fid);
end