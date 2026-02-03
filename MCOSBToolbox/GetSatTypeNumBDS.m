% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [BSatTypeNum,Bsat,BPRN] = GetSatTypeNumBDS(NumRec,Path,list_bds)
    BPRN=linspace(0,0,28);
    for i=1:NumRec
        load([Path '/' list_bds(i).name],'-mat');
        for j=1:28
            if BDSSPRMCB(4,j)~=0
                BPRN(j)=BPRN(j)+1;
            end
        end
        clear BDSSPRMCB;
    end   
    Bsat=find(BPRN==0);
    if isempty(Bsat)
        BSatTypeNum=28;
    else
        % the number of satellites
        BSatTypeNum=28-length(Bsat);
    end
end