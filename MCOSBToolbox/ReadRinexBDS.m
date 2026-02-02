% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function ReadRinexBDS(ipath,opath,STA_opath)
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
        [obs,coor,RCVType]=read_rnx3BDS3(fname);
        % Save observations
        if  ~isempty(obs) && sum(coor) ~= 0
            save(filename1,'obs','coor','-mat');
            STAInformation.coor(i,:)   =coor;
            STAInformation.RCVType{1,i}=RCVType;
        end
        disp(['-------- [',num2str(doy),'] [ ',num2str(i),'/',num2str(len),' ] BDS3OBS file is read !']);
    end
    save(filename2,'STAInformation','-mat');
end
% Read OBS from RINEX3 files
function [obs,coor,RCVType]=read_rnx3BDS3(path)
% INPUT: 
%   path: path of observation files
% OUTPUT: 
%   obs:  struct of observation files
%   coor: station coordinates
%   RCVType: receiver type
BDSFlag = 0;
bdsmark = 0;
fid=fopen(path,'r');
while 1
    line=fgetl(fid);
    if ~ischar(line)
        break
    end
    % the rinex file is null or end
    % RINEX version should be larger than 3.02
    if length(line)>79 && strcmp(line(61:80),'RINEX VERSION / TYPE')
        if strcmp(line(6),'3')==0
            obs=[];coor=[0,0,0]; RCVType='';
            break;
        end
        if str2double(line(9))<3
            obs=[];coor=[0,0,0]; RCVType='';
            break;
        end
    end
    % get receivers coordinates-----------
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
    % get the BDS observable types
    if length(line)>78 && strcmp('C',line(1)) && strcmp('SYS / # / OBS TYPES',line(61:79))
        bdsloc=zeros(1,17);
        obsbds_n=str2double(line(4:6));
        if obsbds_n>39
            for i=1:13
                bdst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                bdst(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                bdst(i+26)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsbds_n-39
                bdst(i+39)={line(4+4*i:6+4*i)};
            end

        elseif obsbds_n>26 && obsbds_n<=39
            for i=1:13
                bdst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                bdst(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsbds_n-26
                bdst(i+26)={line(4+4*i:6+4*i)};
            end

        elseif obsbds_n>13 && obsbds_n<=26
            for i=1:13
                bdst(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsbds_n-13
                bdst(i+13)={line(4+4*i:6+4*i)};
            end

        else
            for i=1:obsbds_n
                bdst(i)={line(4+4*i:6+4*i)};
            end
        end
        % B1I：C2I/C2Q/C2X
        % B1C：C1D/C1P/C1X/C1A
        % B2a：C5D/C5P/C5X
        % B2b：C7D/C7P/C7Z
        % B3I：C6I/C6Q/C6X/C6A
        for i=1:obsbds_n
            if strcmp(bdst(1,i),'C2I')
                bdsloc(1)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C2Q')
                bdsloc(2)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C2X')
                bdsloc(3)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C1D')
                bdsloc(4)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C1P') 
                bdsloc(5)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C1X') 
                bdsloc(6)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C1A')
                bdsloc(7)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C5D')
                bdsloc(8)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C5P')
                bdsloc(9)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C5X')
                bdsloc(10)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C7D')
                bdsloc(11)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C7P')
                bdsloc(12)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C7Z')
                bdsloc(13)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C6I') 
                bdsloc(14)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C6Q') 
                bdsloc(15)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C6X')
                bdsloc(16)=i;
                continue;
            end
            if strcmp(bdst(1,i),'C6A')
                bdsloc(17)=i;
                continue;
            end
        end
        % Note: 
        % the following code types should be observed for BDS3
        % at least five-frequency and C2I/C6I
        % C1P/C2I/C5P/C6I/C7D
        NUML1 = sum(bdsloc(1:3));
        NUML2 = sum(bdsloc(4:7));
        NUML3 = sum(bdsloc(8:10));
        NUML4 = sum(bdsloc(11:13));
        NUML5 = sum(bdsloc(14:17));
        if NUML1>0 && NUML2>0 && NUML3>0 && NUML4>0 && NUML5>0 && bdsloc(1)>0 && bdsloc(14)>0
            if (bdsloc(5)>0&&bdsloc(9)>0) || (bdsloc(6)>0&&bdsloc(10)>0)
                BDSFlag = 1;
            end
        end
        continue;  
    end     
    % start getting observables
    if(length(line)>72&&strcmp(line(61:73),'END OF HEADER')) 
        if BDSFlag == 0
            obs=[];coor=[0,0,0];RCVType='';
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
            if exist('bdsloc','var')
                if strcmp(line(1),'C')
                    bdsmark = 1;
                    sv_C=str2double(line(2:3));
                    if isnan(sv_C) 
                        continue; 
                    end
                    if bdsloc(1)~=0&&length(line)>bdsloc(1)*16    obs.BDSC2I(ep,sv_C)=str2double(line(16*bdsloc(1)-12:16*bdsloc(1)+1)); end
                    if bdsloc(2)~=0&&length(line)>bdsloc(2)*16    obs.BDSC2Q(ep,sv_C)=str2double(line(16*bdsloc(2)-12:16*bdsloc(2)+1)); end
                    if bdsloc(3)~=0&&length(line)>bdsloc(3)*16    obs.BDSC2X(ep,sv_C)=str2double(line(16*bdsloc(3)-12:16*bdsloc(3)+1)); end
                    if bdsloc(4)~=0&&length(line)>bdsloc(4)*16    obs.BDSC1D(ep,sv_C)=str2double(line(16*bdsloc(4)-12:16*bdsloc(4)+1)); end
                    if bdsloc(5)~=0&&length(line)>bdsloc(5)*16    obs.BDSC1P(ep,sv_C)=str2double(line(16*bdsloc(5)-12:16*bdsloc(5)+1)); end
                    if bdsloc(6)~=0&&length(line)>bdsloc(6)*16    obs.BDSC1X(ep,sv_C)=str2double(line(16*bdsloc(6)-12:16*bdsloc(6)+1)); end
                    if bdsloc(7)~=0&&length(line)>bdsloc(7)*16    obs.BDSC1A(ep,sv_C)=str2double(line(16*bdsloc(7)-12:16*bdsloc(7)+1)); end
                    if bdsloc(8)~=0&&length(line)>bdsloc(8)*16    obs.BDSC5D(ep,sv_C)=str2double(line(16*bdsloc(8)-12:16*bdsloc(8)+1)); end
                    if bdsloc(9)~=0&&length(line)>bdsloc(9)*16    obs.BDSC5P(ep,sv_C)=str2double(line(16*bdsloc(9)-12:16*bdsloc(9)+1)); end
                    if bdsloc(10)~=0&&length(line)>bdsloc(10)*16  obs.BDSC5X(ep,sv_C)=str2double(line(16*bdsloc(10)-12:16*bdsloc(10)+1)); end
                    if bdsloc(11)~=0&&length(line)>bdsloc(11)*16  obs.BDSC7D(ep,sv_C)=str2double(line(16*bdsloc(11)-12:16*bdsloc(11)+1)); end
                    if bdsloc(12)~=0&&length(line)>bdsloc(12)*16  obs.BDSC7P(ep,sv_C)=str2double(line(16*bdsloc(12)-12:16*bdsloc(12)+1)); end                  
                    if bdsloc(13)~=0&&length(line)>bdsloc(13)*16  obs.BDSC7Z(ep,sv_C)=str2double(line(16*bdsloc(13)-12:16*bdsloc(13)+1)); end
                    if bdsloc(14)~=0&&length(line)>bdsloc(14)*16  obs.BDSC6I(ep,sv_C)=str2double(line(16*bdsloc(14)-12:16*bdsloc(14)+1)); end
                    if bdsloc(15)~=0&&length(line)>bdsloc(15)*16  obs.BDSC6Q(ep,sv_C)=str2double(line(16*bdsloc(15)-12:16*bdsloc(15)+1)); end
                    if bdsloc(16)~=0&&length(line)>bdsloc(16)*16  obs.BDSC6X(ep,sv_C)=str2double(line(16*bdsloc(16)-12:16*bdsloc(16)+1)); end
                    if bdsloc(17)~=0&&length(line)>bdsloc(17)*16  obs.BDSC6A(ep,sv_C)=str2double(line(16*bdsloc(17)-12:16*bdsloc(17)+1)); end
                end
            end
        end
        if bdsmark == 0
            obs=[];coor=[0,0,0]; RCVType='';           
        end
    end
end
fclose(fid);
end