% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [fid_write] = WriteOSBGAL(GALSatOSBC1C,GALSatOSBC1X,fid_write,index,DOYth,Year)
    SATOSB1=GALSatOSBC1C(index).OSB;
    SATOSB2=GALSatOSBC1X(index).OSB;
    j = 1;
    for i = 1:36     
        C1C = SATOSB1(j).C1C;
        C1X = SATOSB2(j).C1X;
        C6C = SATOSB1(j).C6C;
        C6X = SATOSB2(j).C6X;
        C7Q = SATOSB1(j).C7Q;
        C7X = SATOSB2(j).C7X;
        C8Q = SATOSB1(j).C8Q;
        C8X = SATOSB2(j).C8X;
        C5Q = SATOSB1(j).C5Q;
        C5X = SATOSB2(j).C5X;
        prn = sprintf(' E%02d', i);
        DOY1 = sprintf('%03d', DOYth);
        DOY2 = sprintf('%03d', (DOYth+1));
        str = [num2str(Year) ':' DOY1 ':' '00000 ' num2str(Year) ':' DOY2 ':' '00000 ' 'ns     '];
        if (abs(C1C)>0)
            num = sprintf('%19.4f      %1.4f', C1C, 0.000);
            line = [' OSB  ' 'E000' prn '           C1C       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C1X)>0)
            num = sprintf('%19.4f      %1.4f', C1X, 0.000);
            line = [' OSB  ' 'E000' prn '           C1X       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C6C)>0)
            num = sprintf('%19.4f      %1.4f', C6C, 0.000);
            line = [' OSB  ' 'E000' prn '           C6C       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C6X)>0)
            num = sprintf('%19.4f      %1.4f', C6X, 0.000);
            line = [' OSB  ' 'E000' prn '           C6X       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C7Q)>0)
            num = sprintf('%19.4f      %1.4f', C7Q, 0.000);
            line = [' OSB  ' 'E000' prn '           C7Q       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C7X)>0)
            num = sprintf('%19.4f      %1.4f', C7X, 0.000);
            line = [' OSB  ' 'E000' prn '           C7X       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C8Q)>0)
            num = sprintf('%19.4f      %1.4f', C8Q, 0.000);
            line = [' OSB  ' 'E000' prn '           C8Q       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C8X)>0)
            num = sprintf('%19.4f      %1.4f', C8X, 0.000);
            line = [' OSB  ' 'E000' prn '           C8X       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C5Q)>0)
            num = sprintf('%19.4f      %1.4f', C5Q, 0.000);
            line = [' OSB  ' 'E000' prn '           C5Q       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C5X)>0)
            num = sprintf('%19.4f      %1.4f', C5X, 0.000);
            line = [' OSB  ' 'E000' prn '           C5X       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        clear C1C C1X C6C C6X C7Q C7X C8Q C8X C5Q C5X;
        j = j + 1;
    end
end