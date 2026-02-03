% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [fid_write] = WriteOSBGPS(GPSSatOSBC5Q,GPSSatOSBC5X,fid_write,index,DOYth,Year)
    SATOSB1=GPSSatOSBC5Q(index).OSB;
    SATOSB2=GPSSatOSBC5X(index).OSB;
    j = 1;
    for i = 1:32
        C1C = SATOSB1(j).C1C;
        C1W = SATOSB1(j).C1W;
        C2W = SATOSB1(j).C2W;
        C2L = SATOSB1(j).C2L;
        C2X = SATOSB2(j).C2X;
        C5Q = SATOSB1(j).C5Q;
        C5X = SATOSB2(j).C5X;
        prn = sprintf(' G%02d', i);
        DOY1 = sprintf('%03d', DOYth);
        DOY2 = sprintf('%03d', (DOYth+1));
        str = [num2str(Year) ':' DOY1 ':' '00000 ' num2str(Year) ':' DOY2 ':' '00000 ' 'ns     '];
        if (abs(C1C)>0)
            num = sprintf('%19.4f      %1.4f', C1C, 0.000);
            line = [' OSB  ' 'G000' prn '           C1C       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C1W)>0)
            num = sprintf('%19.4f      %1.4f', C1W, 0.000);
            line = [' OSB  ' 'G000' prn '           C1W       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C2W)>0)
            num = sprintf('%19.4f      %1.4f', C2W, 0.000);
            line = [' OSB  ' 'G000' prn '           C2W       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C2L)>0)
            num = sprintf('%19.4f      %1.4f', C2L, 0.000);
            line = [' OSB  ' 'G000' prn '           C2L       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C2X)>0)
            num = sprintf('%19.4f      %1.4f', C2X, 0.000);
            line = [' OSB  ' 'G000' prn '           C2X       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C5Q)>0)
            num = sprintf('%19.4f      %1.4f', C5Q, 0.000);
            line = [' OSB  ' 'G000' prn '           C5Q       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C5X)>0)
            num = sprintf('%19.4f      %1.4f', C5X, 0.000);
            line = [' OSB  ' 'G000' prn '           C5X       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        clear C1C C1W C2W C2L C2X C5Q C5X;
        j = j + 1;
    end
end