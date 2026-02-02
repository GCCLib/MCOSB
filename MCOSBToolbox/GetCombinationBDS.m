% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [TempBDSSPRMCB] = GetCombinationBDS(Bsat,BDSSPRMCB)
    % Delete zero colum, i.e., without the observed satellites
    TempBDSSPRMCB = [];
    if ~isempty(Bsat)
        for k = 1:28
            if sum(Bsat == k) > 0
               continue;
            else
               TempBDSSPRMCB = [TempBDSSPRMCB BDSSPRMCB(:,k)];
            end
        end
    else
        % full
        TempBDSSPRMCB = BDSSPRMCB;
    end
end