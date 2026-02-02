% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [BDSOSB] = EstimateOSBBDS(NumSat,NumPar,NumRec,NBB,W,NumFre)
    FB1I = 1.561098E9;
    FB3I = 1.26852E9;
    IFa =  FB1I^2/(FB1I^2-FB3I^2);
    IFb = -FB3I^2/(FB1I^2-FB3I^2);
    % Constraint conditions
    % 1: B1I/B3I IF = 0 for each receriver
    C1 = [];
    for j=1:NumRec
        iRecst = NumRec+j;
        iReced = NumRec+NumRec+j;
        CTemp = linspace(0,0,NumPar);
        CTemp(iRecst) = IFa;
        CTemp(iReced) = IFb;
        C1 = [C1;CTemp];
        clear CTemp;
    end
    % 1: B1I/B3I IF = 0 for each satellite
    C2 = C1;
    for j=1:NumSat
        n_Rec = NumRec*NumFre;
        iSatst = n_Rec + NumSat+ j;
        iSated = n_Rec + NumSat+ NumSat + j;
        CTemp = linspace(0,0,NumPar);
        CTemp(iSatst) = IFa;
        CTemp(iSated) = IFb;
        C2 = [C2;CTemp];
        clear CTemp;
    end
    % 3: zero-value for each frequency
    C3 = C2;
    n_Rec = NumRec*NumFre;
    iSatst = n_Rec + 1;
    iSated = n_Rec + NumSat;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,NumSat);
    C3 = [C3;CTemp];
    clear CTemp;

    iSatst = n_Rec +NumSat+ 1;
    iSated = n_Rec +NumSat+ NumSat;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,NumSat);
    C3 = [C3;CTemp];
    clear CTemp;

    iSatst = n_Rec +NumSat+NumSat+ 1;
    iSated = n_Rec +NumSat+NumSat+NumSat;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,NumSat);
    C3 = [C3;CTemp];
    clear CTemp;

    iSatst = n_Rec +NumSat+NumSat+NumSat+ 1;
    iSated = n_Rec +NumSat+NumSat+NumSat+NumSat;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,NumSat);
    C3 = [C3;CTemp];
    clear CTemp;

    iSatst = n_Rec +NumSat+NumSat+NumSat+NumSat+ 1;
    iSated = n_Rec +NumSat+NumSat+NumSat+NumSat+NumSat;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,NumSat);
    C3 = [C3;CTemp];
    clear CTemp;    

    % Estimation
    Wx = zeros(NumRec+NumSat+NumFre,1);
    C = C3;
    NBB=NBB+C'*C;
    W=W+C'*Wx;
    OSB=pinv(NBB)*W;
    % ns
    % RecBDSOSB=OSB(1:NumRecPar)*10^9/299792458;
    BDSOSB=OSB(NumRec*NumFre+1:end)*10^9/299792458;
end