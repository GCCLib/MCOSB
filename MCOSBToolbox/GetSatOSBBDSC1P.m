% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [BDSSatOSB] = GetSatOSBBDSC1P(NumSat,BPRN,BDSOSB)
    % B1C
    count = 0;
    for i=1:28
        if BPRN(1,i)>0
            count = count + 1;
            OSB = BDSOSB(count,1);
        else
            OSB = 0;
        end
        BDSSatOSB(i).C1P = OSB;
    end

    count = 0;
    for i=1:28
        if BPRN(1,i)>0
            count = count + 1;
            OSB = BDSOSB(count+NumSat,1);
        else
            OSB = 0;
        end
        BDSSatOSB(i).C2I = OSB;
    end

    count = 0;
    for i=1:28
        if BPRN(1,i)>0
            count = count + 1;
            OSB = BDSOSB(count+NumSat+NumSat,1);
        else
            OSB = 0;
        end
        BDSSatOSB(i).C6I = OSB;
    end

    count = 0;
    for i=1:28
        if BPRN(1,i)>0
            count = count + 1;
            OSB = BDSOSB(count+NumSat+NumSat+NumSat,1);
        else
            OSB = 0;
        end
        BDSSatOSB(i).C7D = OSB;
    end

    count = 0;
    for i=1:28
        if BPRN(1,i)>0
            count = count + 1;
            OSB = BDSOSB(count+NumSat+NumSat+NumSat+NumSat,1);
        else
            OSB = 0;
        end
        BDSSatOSB(i).C5P = OSB;
    end
end