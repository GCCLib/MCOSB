% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [TempGALSPRMCB] = GetCombinationGAL(Esat,GALSPRMCB)
    % Delete zero colum, i.e., without the observed satellites
    TempGALSPRMCB = [];
    if ~isempty(Esat)
        for k = 1:36
            if sum(Esat == k) > 0
               continue;
            else
               TempGALSPRMCB = [TempGALSPRMCB GALSPRMCB(:,k)];
            end
        end
    else
        % full
        TempGALSPRMCB = GALSPRMCB;
    end
end