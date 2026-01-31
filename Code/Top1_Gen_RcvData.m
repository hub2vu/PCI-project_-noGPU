%  2019-04-01
%     Written by Sua Bae
%     Ver 1.0: Gen Simulation Data (mat and bin files): BW 100%, 7 point targets
%  2025-09-04
%     Modified by Sua Bae
% 
clear; close all;
addpath(genpath('src'));

% transducer
load('stTrans_L7-4.mat');
stTrans.nBW_6dB = 100; % [%] -6dB bandwidth of transducer

% data type
exp_idx = 1;        % 1:simulation, 2:phantom, 3:invivo
phantom_idx = 3;
acqType_idx = 2;    % 1:CF, 2:PW
ang_deg = 5;       % 0 or 10 (deg of max PW ang)
             
% directory
stD = clDir(exp_idx, phantom_idx, acqType_idx, ang_deg, '..\Data\');
CreateFolders(stD);

% Tx
stTx.sType = 'PW';
stTx.nNum = 21;
stTx.aPwAngle_deg = linspace(-ang_deg,ang_deg,stTx.nNum);
stTx.nFc = stTrans.nFc;
stTx.nTxPulseCycle = 1;

% Scatterers
if phantom_idx < 7
    stScat.mScatXYZPos = [0e-3,0, phantom_idx*10e-3];
elseif phantom_idx == 7
    stScat.mScatXYZPos = [0e-3,0,10e-3;
                          0e-3,0,20e-3;
                          0e-3,0,30e-3;
                          0e-3,0,40e-3;
                          0e-3,0,50e-3;
                          0e-3,0,60e-3];
end
stScat.aScatMag = 1e26*ones(size(stScat.mScatXYZPos,1),1);
    
bPlot = 1;

% Call function:: GenRFdata
stRcvData = GenRFdata(stTrans, stTx, stScat, bPlot);

hFig = figure; hAx = axes;
% save mRcvData
for txidx = 1:stRcvData.stTx.nNum
    mRcvData = stRcvData.vRcvData(:,:,txidx);
    save([stD.sDataDir 'RcvData\' 'mRcvData_tx_' num2str(txidx) '.mat'],'mRcvData');
    imagesc(hAx, mRcvData); title(['mRcvData-tx-' num2str(txidx) ', angle=' num2str(stTx.aPwAngle_deg(txidx)) ' deg']); pause(0.1);
end


% save stRFInfo
stRFInfo.sVersion       = stRcvData.sVersion; % '2.1' stRFInfo.nRcvChannel included, TxAptrGen.m used
stRFInfo.stTrans        = stTrans;
stRFInfo.nRcvChannel    = 128;  
stRFInfo.nSoundSpeed    = stRcvData.nSoundSpeed; % (m/s)
stRFInfo.nSample        = size(mRcvData,1);
stRFInfo.nFs            = stRcvData.nFs;
stRFInfo.nFc            = stRcvData.nFc;
stRFInfo.stTx           = stRcvData.stTx;
stRFInfo.stScat         = stRcvData.stScat;
stRFInfo.sStartAcqTime  = stRcvData.sStartAcqTime;
save([stD.sDataDir 'RcvData\' 'stRFInfo.mat'],'stRFInfo');
