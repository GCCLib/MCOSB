% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [GRecTypeNum] = GetRecTypeNumC5X(Path)
    list_obs=dir([Path,'/*.mat']);
    len=length(list_obs);
    GPSC1C = 0;
    GPSC1W = 0;
    GPSC2W = 0;
    GPSC2X = 0;
    GPSC5X = 0;
    for i=1:len
        load([Path,'/',list_obs(i).name],'-mat');
        if sum(GPSSPRMCB(1,:))~=0
            GPSC1W = GPSC1W + 1;
            GPSC2W = GPSC2W + 1;
        end        
        if sum(GPSSPRMCB(3,:))~=0
            GPSC5X = GPSC5X + 1;
        end        
        if sum(GPSSPRMCB(4,:))~=0
            GPSC1C = GPSC1C + 1;
        end
        if sum(GPSSPRMCB(6,:))~=0
            GPSC2X = GPSC2X + 1;
        end
        clear GPSSPRMCB;
    end
   GRecTypeNum = [GPSC1C GPSC1W GPSC2W GPSC2X GPSC5X];
end