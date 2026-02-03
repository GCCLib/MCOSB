% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [GRecTypeNum] = GetRecTypeNumC5Q(Path)
    list_obs=dir([Path,'/*.mat']); 
    len=length(list_obs);
    GPSC1C = 0;
    GPSC1W = 0;
    GPSC2W = 0;
    GPSC2L = 0;
    GPSC5Q = 0;   
    for i = 1:len
        load([Path,'/',list_obs(i).name],'-mat');
        if sum(GPSSPRMCB(1,:))~=0
            GPSC1W = GPSC1W + 1;
            GPSC2W = GPSC2W + 1;
        end
        if sum(GPSSPRMCB(2,:))~=0
            GPSC5Q = GPSC5Q + 1;
        end 
        if sum(GPSSPRMCB(4,:))~=0
            GPSC1C = GPSC1C + 1;
        end
        if sum(GPSSPRMCB(5,:))~=0
            GPSC2L = GPSC2L + 1;
        end
        if sum(GPSSPRMCB(5,:))==0
            debug = 1;
        end
        clear GPSSPRMCB;
    end
    GRecTypeNum = [GPSC1C GPSC1W GPSC2W GPSC2L GPSC5Q];
end