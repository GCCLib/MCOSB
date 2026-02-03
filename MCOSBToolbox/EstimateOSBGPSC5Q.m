% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [GPSOSB,SatC1C,SatC1W,SatC2W,SatC2L] = EstimateOSBGPSC5Q(GSatTypeNum,RecTypeNum,NumPar,NumRecPar,NumRec,NBB,W)
    %Number of code type
    NumType = 5;
    FL1  = 1.57542E9;
    FL2  = 1.22760E9;
    IFa =  FL1^2/(FL1^2-FL2^2);
    IFb = -FL2^2/(FL1^2-FL2^2);
    % number of satellite OSB
    SatC1C=GSatTypeNum(1);
    SatC1W=GSatTypeNum(2);
    SatC2W=GSatTypeNum(3);
    SatC2L=GSatTypeNum(4);
    SatC5Q=GSatTypeNum(5);
    % number of receiver OSB
    RecC1C=RecTypeNum(1);
    RecC1W=RecTypeNum(2);
    % Constraint conditions
    % 1: C1W/C2W IF = 0 for each receriver
    C1 = [];
    for j = 1:RecC1W
        iRecst = RecC1C + j;
        iReced = RecC1C + RecC1W + j;
        CTemp = linspace(0,0,NumPar);
        CTemp(iRecst) = IFa;
        CTemp(iReced) = IFb;
        C1 = [C1;CTemp];
        clear CTemp;
    end
    % 2: C1W/C2W IF = 0 for each satellite
    C2 = C1;
    for j=1:SatC1W
        iSatst = NumRecPar + SatC1C + j;
        iSated = NumRecPar + SatC1C + SatC1W + j;
        CTemp = linspace(0,0,NumPar);
        CTemp(iSatst) = IFa;
        CTemp(iSated) = IFb;
        C2 = [C2;CTemp];
        clear CTemp;
    end
    % 3: zero-value for each frequency
    % C1C
    C3 = C2;
    iSatst = NumRecPar + 1;
    iSated = NumRecPar + SatC1C;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,SatC1C);
    C3 = [C3;CTemp];
    clear CTemp;
    % C1W
    iSatst = NumRecPar +SatC1C+ 1;
    iSated = NumRecPar +SatC1C + SatC1W;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,SatC1W);
    C3 = [C3;CTemp];
    clear CTemp;
    % C2W
    iSatst = NumRecPar +SatC1C + SatC1W + 1;
    iSated = NumRecPar +SatC1C + SatC1W + SatC2W;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,SatC2W);
    C3 = [C3;CTemp];
    clear CTemp;
    % C2L
    iSatst = NumRecPar +SatC1C + SatC1W + SatC2W + 1;
    iSated = NumRecPar +SatC1C + SatC1W + SatC2W + SatC2L;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,SatC2L);
    C3 = [C3;CTemp];
    clear CTemp;
    % C5Q
    iSatst = NumRecPar +SatC1C + SatC1W + SatC2W + SatC2L + 1;
    iSated = NumRecPar +SatC1C + SatC1W + SatC2W + SatC2L + SatC5Q;
    CTemp = linspace(0,0,NumPar);
    CTemp(iSatst:iSated) = ones(1,SatC5Q);
    C3 = [C3;CTemp];
    clear CTemp;
    % Estimation
    Wx = zeros(NumRec+SatC1W+NumType,1);
    C = C3;
    NBB=NBB+C'*C;
    W=W+C'*Wx;
    OSB=pinv(NBB)*W;
    % ns
    % RecGPSOSB=OSB(1:NumRecPar)*10^9/299792458;
    GPSOSB=OSB(NumRecPar+1:end)*10^9/299792458;
end