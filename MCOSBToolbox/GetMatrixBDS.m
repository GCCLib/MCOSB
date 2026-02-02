% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [M,l,P]=GetMatrixBDS(SPRMCB,NumRec,NumSat,ith,NumPar,NumFre)
    % Design matrix of parameters
    M=[]; l=[]; P=[];
    % Constants
    CLIGHT = 299792458;
    FB1C = 1.57542E9;
    FB1I = 1.561098E9;
    FB3I = 1.26852E9;
    FB2B = 1.20714E9;
    FB2A = 1.17645E9;
    % 1
    WLB1C = CLIGHT/FB1C; 
    % 2
    WLB1I = CLIGHT/FB1I;
    % 3
    WLB3I = CLIGHT/FB3I; 
    % 4
    WLB2B = CLIGHT/FB2B;
    % 5
    WLB2A = CLIGHT/FB2A;
    A32 = SQR(WLB3I) - SQR(WLB1I);
    A13 = SQR(WLB1C) - SQR(WLB3I);
    A21 = SQR(WLB1I) - SQR(WLB1C);
    A43 = SQR(WLB2B) - SQR(WLB3I);
    A53 = SQR(WLB2A) - SQR(WLB3I);
    A24 = SQR(WLB1I) - SQR(WLB2B);
    A25 = SQR(WLB1I) - SQR(WLB2A);
    % B1I/B3I
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
        M_col(NumRec+ith)=1;
        M_col(NumRec+NumRec+ith)=-1;
        % Satellite
        % i denotes which satellite
        M_col(NumRec*NumFre+NumSat+i)=1;
        M_col(NumRec*NumFre+NumSat+NumSat+i)=-1;
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
        %3 combinations
        for k=1:3
            M_col=zeros(1,NumPar);
            if k == 1
                % B1CB1IB3I
                M_col(ith)=A32;
                M_col(NumRec+ith)=A13;
                M_col(NumRec+NumRec+ith)=A21;
                M_col(NumRec*NumFre+i)=A32;
                M_col(NumRec*NumFre+NumSat+i)=A13;
                M_col(NumRec*NumFre+NumSat+NumSat+i)=A21;
                M=[M;M_col]; 
                l=[l;SPRMCB(k,i)];
                P=[P;1/SPRMCB(4+k,i)];
            end
            if k == 2
                % B1IB3IB2B
                M_col(NumRec+ith)=A43;
                M_col(NumRec+NumRec+ith)=A24;
                M_col(NumRec+NumRec+NumRec+ith)=A32;
                M_col(NumRec*NumFre+NumSat+i)=A43;
                M_col(NumRec*NumFre+NumSat+NumSat+i)=A24;
                M_col(NumRec*NumFre+NumSat+NumSat+NumSat+i)=A32;
                M=[M;M_col]; 
                l=[l;SPRMCB(k,i)];
                P=[P;1/SPRMCB(4+k,i)];
            end
            if k == 3
                % B1IB3IB2A
                M_col(NumRec+ith)=A53;
                M_col(NumRec+NumRec+ith)=A25;
                M_col(NumRec+NumRec+NumRec+NumRec+ith)=A32;
                M_col(NumRec*NumFre+NumSat+i)=A53;
                M_col(NumRec*NumFre+NumSat+NumSat+i)=A25;
                M_col(NumRec*NumFre+NumSat+NumSat+NumSat+NumSat+i)=A32;
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