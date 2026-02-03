% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [fid_write] = WriteOSBBDS(BDSSatOSBC1P,BDSSatOSBC1X,fid_write,index,DOYth,Year)
    SATOSB1=BDSSatOSBC1P(index).OSB;
    SATOSB2=BDSSatOSBC1X(index).OSB;
    j = 1;
    for i = 1:46
        if i < 19
            continue;
        end
        C1P = SATOSB1(j).C1P;
        C1X = SATOSB2(j).C1X;
        C2I = SATOSB1(j).C2I;
        C6I = SATOSB1(j).C6I;
        C7D = SATOSB1(j).C7D;
        C7Z = SATOSB2(j).C7Z;
        C5P = SATOSB1(j).C5P;
        C5X = SATOSB2(j).C5X;
        prn = sprintf(' C%02d', i);
        DOY1 = sprintf('%03d', DOYth);
        DOY2 = sprintf('%03d', (DOYth+1));
        str = [num2str(Year) ':' DOY1 ':' '00000 ' num2str(Year) ':' DOY2 ':' '00000 ' 'ns     '];
        if (abs(C1P)>0)
            num = sprintf('%19.4f      %1.4f', C1P, 0.000);
            line = [' OSB  ' 'C000' prn '           C1P       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C1X)>0)
            num = sprintf('%19.4f      %1.4f', C1X, 0.000);
            line = [' OSB  ' 'C000' prn '           C1X       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C2I)>0)
            num = sprintf('%19.4f      %1.4f', C2I, 0.000);
            line = [' OSB  ' 'C000' prn '           C2I       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C6I)>0)
            num = sprintf('%19.4f      %1.4f', C6I, 0.000);
            line = [' OSB  ' 'C000' prn '           C6I       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C7D)>0)
            num = sprintf('%19.4f      %1.4f', C7D, 0.000);
            line = [' OSB  ' 'C000' prn '           C7D       ' str num];
            fprintf(fid_write, '%s\n', line); 
        end
        if (abs(C7Z)>0)
            num = sprintf('%19.4f      %1.4f', C7Z, 0.000);
            line = [' OSB  ' 'C000' prn '           C7Z       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C5P)>0)
            num = sprintf('%19.4f      %1.4f', C5P, 0.000);
            line = [' OSB  ' 'C000' prn '           C5P       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        if (abs(C5X)>0)
            num = sprintf('%19.4f      %1.4f', C5X, 0.000);
            line = [' OSB  ' 'C000' prn '           C5X       ' str num];
            fprintf(fid_write, '%s\n', line);
        end
        clear C1P C1X C2I C6I C7D C7Z C5P C5X;
        j = j + 1;
    end
    line = '-BIAS/SOLUTION';
    fprintf(fid_write, '%s\n', line);
    line = '%=ENDBIA';
    fprintf(fid_write, '%s\n', line);                                                                   
end