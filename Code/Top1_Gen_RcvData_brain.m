% Top1_Gen_RcvData_brain.m
% Brain atlas annotation -> scatterer cloud -> RF data generation (Field II)

clear; close all; clc;

addpath(genpath('src'));
addpath('simulation');

% transducer
load('stTrans_L7-4.mat');
stTrans.nBW_6dB = 100; % [%]

% data type
exp_idx = 1;        % 1:simulation
phantom_idx = 3;    % clDir/폴더 구조 깨지기 싫어서 기존 값 유지
acqType_idx = 2;    % 1:CF, 2:PW
ang_deg = 4;        % Top3와 반드시 맞추기 (여기 4로 통일)

% directory (기존 Data와 섞이지 않게 별도 폴더 추천)
stD = clDir(exp_idx, phantom_idx, acqType_idx, ang_deg, '..\Data_brain\');
CreateFolders(stD);

% Tx (PW)
stTx.sType = 'PW';
stTx.nNum = 21;
stTx.aPwAngle_deg = linspace(-ang_deg, ang_deg, stTx.nNum);
stTx.nFc = stTrans.nFc;
stTx.nTxPulseCycle = 1;

% ===== Brain Scatterers =====
NRRD_PATH = fullfile('simulation', 'annotation_25.nrrd');
N_SCAT  = 200000;  % 5e4~3e5
SLAB_MM = 0.4;     % 2D 단면 두께
Z0_MM   = 0.6;     % standoff

% amplitude 튜닝(너무 하얗게 뜨면 amp_scale 줄이기)
AMP_SCALE = 1e26;
AMP_SIGMA = 0.6;

stScat = make_brain_scat_from_annotation(NRRD_PATH, N_SCAT, SLAB_MM, Z0_MM, AMP_SCALE, AMP_SIGMA);

bPlot = 0;

% Call function:: GenRFdata
stRcvData = GenRFdata(stTrans, stTx, stScat, bPlot);

% save mRcvData (각 tx별 저장)
hFig = figure('Color','w'); hAx = axes(hFig);
for txidx = 1:stRcvData.stTx.nNum
    mRcvData = stRcvData.vRcvData(:,:,txidx);
    save([stD.sDataDir 'RcvData\' 'mRcvData_tx_' num2str(txidx) '.mat'],'mRcvData');
    imagesc(hAx, mRcvData);
    title(['mRcvData tx=' num2str(txidx) ', angle=' num2str(stTx.aPwAngle_deg(txidx)) ' deg']);
    drawnow;
end

% save stRFInfo
stRFInfo.sVersion       = stRcvData.sVersion; % '2.1' etc
stRFInfo.stTrans        = stTrans;
stRFInfo.nRcvChannel    = 128;
stRFInfo.nSoundSpeed    = stRcvData.nSoundSpeed;
stRFInfo.nSample        = size(mRcvData,1);
stRFInfo.nFs            = stRcvData.nFs;
stRFInfo.nFc            = stRcvData.nFc;
stRFInfo.stTx           = stRcvData.stTx;
stRFInfo.stScat         = stRcvData.stScat;
stRFInfo.sStartAcqTime  = stRcvData.sStartAcqTime;

save([stD.sDataDir 'RcvData\' 'stRFInfo.mat'],'stRFInfo');

disp("✅ Brain RF generation done.");
disp(['Saved to: ' stD.sDataDir]);
