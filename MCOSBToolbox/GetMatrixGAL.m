% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [M,l,P]=GetMatrixGAL(SPRMCB,NumRec,NumSat,ith,NumPar,NumFre)
    % Design matrix of parameters
    M=[]; l=[]; P=[];
    % Constants
    CLIGHT = 299792458;   
    FE1 = 1.57542E9;
    FE6 = 1.27875E9;
    FE5B = 1.20714E9;
    FE5 = 1.191795E9;
    FE5A = 1.17645E9;
    % 1
    WLE1 = CLIGHT/FE1; 
    % 2
    WLE6 = CLIGHT/FE6;
    % 3
    WLE5B = CLIGHT/FE5B; 
    % 4
    WLE5 = CLIGHT/FE5;
    % 5
    WLE5A = CLIGHT/FE5A;    
    A21 = SQR(WLE6) - SQR(WLE1);
    A52 = SQR(WLE5A) - SQR(WLE6);
    A15 = SQR(WLE1) - SQR(WLE5A);
    A31 = SQR(WLE5B) - SQR(WLE1);
    A53 = SQR(WLE5A) - SQR(WLE5B);
    A54 = SQR(WLE5A) - SQR(WLE5);
    A41 = SQR(WLE5) - SQR(WLE1);
    % E1/E5a
    for i=1:NumSat
        if SPRMCB(4,i) == 0
            l=[l;0]; P=[P;0.00000000001];
            M_col=zeros(1,NumPar);
            M=[M;M_col];            
            continue;
        end
        M_col=zeros(1,NumPar);
        % Receiver
        % ith denotes which receiver
        M_col(ith)=1;
        M_col(NumRec+NumRec+NumRec+NumRec+ith)=-1;
        % Satellite
        % i denotes which satellite
        M_col(NumRec*NumFre+i)=1;
        M_col(NumRec*NumFre+NumSat+NumSat+NumSat+NumSat+i)=-1;
        M=[M;M_col];
        l=[l;SPRMCB(4,i)];
        P=[P;1/SPRMCB(8,i)];
    end

    for i=1:NumSat
        if SPRMCB(1,i) == 0 || SPRMCB(2,i) == 0 || SPRMCB(3,i) == 0
            continue;
        end
        if SPRMCB(4,i) == 0
            l=[l;zeros(3,1)];
            P=[P;zeros(3,1)];
            M_col=zeros(3,NumPar);
            M=[M;M_col];
            continue;
        end
        % 3 combinations
        for k=1:3
            M_col=zeros(1,NumPar);
            if k == 1
                M_col(ith)=A52;
                M_col(NumRec+ith)=A15;
                M_col(NumRec+NumRec+NumRec+NumRec+ith)=A21;
                M_col(NumRec*NumFre+i)=A52;
                M_col(NumRec*NumFre+NumSat+i)=A15;
                M_col(NumRec*NumFre+NumSat+NumSat+NumSat+NumSat+i)=A21;
                M=[M;M_col]; 
                l=[l;SPRMCB(k,i)];
                P=[P;1/SPRMCB(4+k,i)];          
            end
            if k == 2
                M_col(ith)=A53;
                M_col(NumRec+NumRec+ith)=A15;
                M_col(NumRec+NumRec+NumRec+NumRec+ith)=A31;
                M_col(NumRec*NumFre+i)=A53;
                M_col(NumRec*NumFre+NumSat+NumSat+i)=A15;
                M_col(NumRec*NumFre+NumSat+NumSat+NumSat+NumSat+i)=A31;
                M=[M;M_col]; 
                l=[l;SPRMCB(k,i)];
                P=[P;1/SPRMCB(4+k,i)];          
            end
            if k == 3
                M_col(ith)=A54;
                M_col(NumRec+NumRec+NumRec+ith)=A15;
                M_col(NumRec+NumRec+NumRec+NumRec+ith)=A41;
                M_col(NumRec*NumFre+i)=A54;
                M_col(NumRec*NumFre+NumSat+NumSat+NumSat+i)=A15;
                M_col(NumRec*NumFre+NumSat+NumSat+NumSat+NumSat+i)=A41;
                M=[M;M_col]; 
                l=[l;SPRMCB(k,i)];
                P=[P;1/SPRMCB(4+k,i)];          
            end     
        end 
    end
    P = diag(P);
end
% SQR
function [outputArg1] = SQR(inputArg1)
    outputArg1 = inputArg1 * inputArg1;
end