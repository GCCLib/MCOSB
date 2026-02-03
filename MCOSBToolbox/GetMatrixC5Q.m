% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [M,l,P]=GetMatrixC5Q(GSatTypeNum,GRecTypeNum,C1WC2W,C1WC2WC5Q,C1WC1C,C2WC2L,ith,NumPar,NumRecPar)
    % Design matrix of parameters
    M=[]; l=[]; P=[];
    % Number of satellite OSB
    SatC1C=GSatTypeNum(1);
    SatC1W=GSatTypeNum(2);
    SatC2W=GSatTypeNum(3);
    SatC2L=GSatTypeNum(4);
    SatC5Q=GSatTypeNum(5);
    % Number of receiver OSB
    RecC1C=GRecTypeNum(1);
    RecC1W=GRecTypeNum(2);
    RecC2W=GRecTypeNum(3);
    RecC2L=GRecTypeNum(4);
    % Constants
    CLIGHT = 299792458;   
    FL1  = 1.57542E9;
    FL2  = 1.22760E9;
    FL5 = 1.17645E9;
    WL1 = CLIGHT/FL1;
    WL2 = CLIGHT/FL2;
    WL5 = CLIGHT/FL5;
    A32 = SQR(WL5) - SQR(WL2);
    A13 = SQR(WL1) - SQR(WL5);
    A21 = SQR(WL2) - SQR(WL1);   
    % C1C/C1W/C2W/C2L/C5Q
    % C1W/C2W
    for i = 1:SatC1W
        % Quality control
        if C1WC2W(1,i)==0 || C1WC2W(7,i) > 0.6
            l=[l;0]; P=[P;0.00000000001];
            M_col=zeros(1,NumPar);
            M=[M;M_col];
            continue;
        end
        M_col=zeros(1,NumPar);
        % Receiver
        % ith denotes which receiver
        M_col(RecC1C+ith) = 1;
        M_col(RecC1C+RecC1W+ith) = -1;
        % Satellite
        % i denotes which satellite
        M_col(NumRecPar+SatC1C+i) = 1;
        M_col(NumRecPar+SatC1C+SatC1W+i) = -1;
        M=[M;M_col];
        l=[l;C1WC2W(1,i)];
        P=[P;1/C1WC2W(7,i)];
    end
    % C1W/C2W/C5Q
    for i=1:SatC5Q
        if C1WC2WC5Q(2,i) == 0
            l=[l;0]; P=[P;0.00000000001];
            M_col=zeros(1,NumPar);
            M=[M;M_col];
            continue;
        end
        M_col=zeros(1,NumPar);
        M_col(RecC1C+ith)=A32;
        M_col(RecC1C+RecC1W+ith)=A13;
        M_col(RecC1C+RecC1W+RecC2W+RecC2L+ith)=A21;
        M_col(NumRecPar+SatC1C+i)=A32;
        M_col(NumRecPar+SatC1C+SatC1W+i)=A13;
        M_col(NumRecPar+SatC1C+SatC1W+SatC2W+SatC2L+i)=A21;
        M=[M;M_col]; 
        l=[l;C1WC2WC5Q(2,i)];
        P=[P;1/C1WC2WC5Q(8,i)];
    end
    % C1W/C1C
    for i=1:SatC1C
        if C1WC1C(4,i) == 0
            l=[l;0]; P=[P;0.00000000001];
            M_col=zeros(1,NumPar);
            M=[M;M_col];
            continue;
        end
        M_col=zeros(1,NumPar);
        M_col(ith)=-1;
        M_col(RecC1C+ith)=1;
        M_col(NumRecPar+i)=-1;
        M_col(NumRecPar+SatC1C+i)=1;
        M=[M;M_col]; 
        l=[l;C1WC1C(4,i)];
        P=[P;1/C1WC1C(10,i)];
    end
    % C2W/C2L
    for i=1:SatC2L
        if C2WC2L(5,i) == 0
            l=[l;0]; P=[P;0.00000000001];
            M_col=zeros(1,NumPar);
            M=[M;M_col];
            continue;
        end
        M_col=zeros(1,NumPar);
        M_col(RecC1C+RecC1W+ith)=1;
        M_col(RecC1C+RecC1W+RecC2W+ith)=-1;
        M_col(NumRecPar+SatC1C+SatC1W+i)=1;
        M_col(NumRecPar+SatC1C+SatC1W+SatC2W+i)=-1;
        M=[M;M_col]; 
        l=[l;C2WC2L(5,i)];
        P=[P;1/C2WC2L(11,i)];
    end
    P = diag(P);
end
% SQR
function [outputArg1] = SQR(inputArg1)
    outputArg1 = inputArg1 * inputArg1;
end