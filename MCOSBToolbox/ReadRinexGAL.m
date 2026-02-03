% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function ReadRinexGAL(ipath,opath,STA_opath)
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
        [obs,coor,RCVType]=read_rnx3GAL(fname); 
        % Save observations
        if  ~isempty(obs) && sum(coor) ~= 0
            save(filename1,'obs','coor','-mat');
            STAInformation.coor(i,:)   =coor;
            STAInformation.RCVType{1,i}=RCVType;
        end
        disp(['-------- [',num2str(doy),'] [ ',num2str(i),'/',num2str(len),' ] GALOBS file is read !']);
    end
    save(filename2,'STAInformation','-mat');
end
% Read OBS from RINEX3 files
function [obs,coor,RCVType]=read_rnx3GAL(path)
% INPUT: 
%   path: path of observation files
% OUTPUT: 
%   obs:  struct of observation files
%   coor: station coordinates
%   RCVType: receiver type
GALFlag = 0;
galmark = 0;
fid=fopen(path,'r');
while 1
    line=fgetl(fid);
    if ~ischar(line), break, end  
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
    % get the GAL observable types
    if length(line)>78 && strcmp('E',line(1)) && strcmp('SYS / # / OBS TYPES',line(61:79))
        GALloc=zeros(1,10);
        obsbds_n=str2double(line(4:6));
        if obsbds_n>39
            for i=1:13
                GALt(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                GALt(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                GALt(i+26)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsbds_n-39
                GALt(i+39)={line(4+4*i:6+4*i)};
            end

        elseif obsbds_n>26 && obsbds_n<=39
            for i=1:13
                GALt(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:13
                GALt(i+13)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsbds_n-26
                GALt(i+26)={line(4+4*i:6+4*i)};
            end

        elseif obsbds_n>13 && obsbds_n<=26
            for i=1:13
                GALt(i)={line(4+4*i:6+4*i)};
            end
            line=fgetl(fid);
            for i=1:obsbds_n-13
                GALt(i+13)={line(4+4*i:6+4*i)};
            end

        else
            for i=1:obsbds_n
                GALt(i)={line(4+4*i:6+4*i)};
            end
        end
        % E1：C1C/C1X
        % E5a：C5Q/C5X
        % E5b：C7Q/C7X
        % E5：C8Q/C8X
        % E6：C6C/C6X
        for i=1:obsbds_n
            if strcmp(GALt(1,i),'C1C')
                GALloc(1)=i;
                continue;
            end
            if strcmp(GALt(1,i),'C1X')
                GALloc(2)=i;
                continue;
            end
            if strcmp(GALt(1,i),'C5Q')
                GALloc(3)=i;
                continue;
            end
            if strcmp(GALt(1,i),'C5X')
                GALloc(4)=i;
                continue;
            end
            if strcmp(GALt(1,i),'C7Q') 
                GALloc(5)=i;
                continue;
            end
            if strcmp(GALt(1,i),'C7X') 
                GALloc(6)=i;
                continue;
            end
            if strcmp(GALt(1,i),'C8Q')
                GALloc(7)=i;
                continue;
            end
            if strcmp(GALt(1,i),'C8X')
                GALloc(8)=i;
                continue;
            end
            if strcmp(GALt(1,i),'C6C')
                GALloc(9)=i;
                continue;
            end
            if strcmp(GALt(1,i),'C6X')
                GALloc(10)=i;
                continue;
            end
        end
        % Note: 
        % the following code types should be observed for BDS3
        % at least five-frequency
        % C1C/C5Q
        % C1X/C5X
        NUML1 = sum(GALloc(1:2));
        NUML2 = sum(GALloc(3:4));
        NUML3 = sum(GALloc(5:6));
        NUML4 = sum(GALloc(7:8));
        NUML5 = sum(GALloc(9:10));
        if NUML1>0 && NUML2>0 && NUML3>0 && NUML4>0 && NUML5>0
            if (GALloc(2)>0&&GALloc(4)>0) || (GALloc(1)>0&&GALloc(3)>0)
                GALFlag = 1;
            end
        end
        continue;  
    end    
    % start getting observables
    if(length(line)>72&&strcmp(line(61:73),'END OF HEADER')) 
        if GALFlag == 0
            obs=[];coor=[0,0,0]; RCVType='';
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
            if exist('GALloc','var')
                if strcmp(line(1),'E')
                    galmark = 1;
                    sv_E=str2double(line(2:3));
                    if isnan(sv_E) 
                        continue; 
                    end
                    if GALloc(1)~=0&&length(line)>GALloc(1)*16    obs.GALC1C(ep,sv_E)=str2double(line(16*GALloc(1)-12:16*GALloc(1)+1)); end
                    if GALloc(2)~=0&&length(line)>GALloc(2)*16    obs.GALC1X(ep,sv_E)=str2double(line(16*GALloc(2)-12:16*GALloc(2)+1)); end
                    if GALloc(3)~=0&&length(line)>GALloc(3)*16    obs.GALC5Q(ep,sv_E)=str2double(line(16*GALloc(3)-12:16*GALloc(3)+1)); end
                    if GALloc(4)~=0&&length(line)>GALloc(4)*16    obs.GALC5X(ep,sv_E)=str2double(line(16*GALloc(4)-12:16*GALloc(4)+1)); end
                    if GALloc(5)~=0&&length(line)>GALloc(5)*16    obs.GALC7Q(ep,sv_E)=str2double(line(16*GALloc(5)-12:16*GALloc(5)+1)); end
                    if GALloc(6)~=0&&length(line)>GALloc(6)*16    obs.GALC7X(ep,sv_E)=str2double(line(16*GALloc(6)-12:16*GALloc(6)+1)); end
                    if GALloc(7)~=0&&length(line)>GALloc(7)*16    obs.GALC8Q(ep,sv_E)=str2double(line(16*GALloc(7)-12:16*GALloc(7)+1)); end
                    if GALloc(8)~=0&&length(line)>GALloc(8)*16    obs.GALC8X(ep,sv_E)=str2double(line(16*GALloc(8)-12:16*GALloc(8)+1)); end
                    if GALloc(9)~=0&&length(line)>GALloc(9)*16    obs.GALC6C(ep,sv_E)=str2double(line(16*GALloc(9)-12:16*GALloc(9)+1)); end
                    if GALloc(10)~=0&&length(line)>GALloc(10)*16  obs.GALC6X(ep,sv_E)=str2double(line(16*GALloc(10)-12:16*GALloc(10)+1)); end
                end
            end
        end
        if galmark == 0
            obs=[];coor=[0,0,0]; RCVType='';           
        end        
    end
end
fclose(fid);
end