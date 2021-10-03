clc;
clear;
load('AngularResolution_12deg_sum_LOS_fewMPCs.mat');
load('UWB_AntennaDiagramm_SGA.mat'); 
RXpredoff= 5;
TXpredoff= RXpredoff; % resolution of predictions of AoD
RXrealoff= 12;% angles of rotation of RX for measurements.
TXrealoff= 12;% angles of rotation of TX for measurements.
RXrealf= 0; % First angle of measurement at RX
RXreall= 359; %last angle of measurement at RX
TXrealf= 0; % First angle of measurement at TX
TXreall= 359; % Last angle of measurement at TX
TXpredf= 0; % First angle of prediction at TX
TXpredl= 359; % Last angle of prediction at TX
RXpredf= 0; % First angle of prediction at RX
RXpredl= 359; % Last angle of prediction at RX
phiRXSteps = fix((RXreall-RXrealf)/RXrealoff)+1; %Number of angles
phiTXSteps = fix((TXreall-TXrealf)/TXrealoff)+1;
phiTXSteps2 = fix((TXpredl-TXpredf)/TXpredoff)+1;
phiRXSteps2 = fix((RXpredl-RXpredf)/RXpredoff)+1;
PAS=reshape(PAPmatrix',phiTXSteps*phiRXSteps,1);
Amat= zeros(phiTXSteps*phiRXSteps, phiTXSteps2*phiRXSteps2);

%% Initialising the sensing matrix used for reconstruction (Amat)
for phiTXpred=TXpredf:TXpredoff:TXpredl
    for phiRXpred= RXpredf:RXpredoff:RXpredl
        for phiTXreal= TXrealf:TXrealoff:TXreall
            for phiRXreal= RXrealf:RXrealoff:RXreall
                Amat(phiRXSteps*floor((phiTXreal-TXrealf)/TXrealoff)+floor((phiRXreal-RXrealf)/RXrealoff)+1, phiRXSteps2*floor((phiTXpred-TXpredf)/TXpredoff)+floor((phiRXpred-RXpredf)/RXpredoff)+1)= pattern_3D.Gain_dBi(181,mod(2*(phiTXpred-phiTXreal)+360, 720)+1)+pattern_3D.Gain_dBi(181,mod(2*(phiRXpred-phiRXreal)+360, 720)+1);
            end
        end
    end
end
Amat=convertdBtoStandard(Amat);

%% Solving y=Ax
sigma=0.00001;
epsilon=1e-6;
b=new_SBL_algo(sigma,Amat,PAS,epsilon); % Reconstruction using l1_ls_nonneg package.
save('CompSens12_5deg_sum.mat', 'PAS', 'b');
