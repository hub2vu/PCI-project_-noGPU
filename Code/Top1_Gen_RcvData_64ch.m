%  2019-04-01
%     Written by Sua Bae
%     Ver 1.0: Gen Simulation Data (mat and bin files): BW 100%, 7 point targets
%  2025-09-04
%     Modified by Sua Bae
%  2025-01-31
%     Modified for 64 channel (skip every other channel from 128ch)
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
ang_deg = 4;        % 0 or 10 (deg of max PW ang)
             
% directory - 64ch 데이터는 Data_64ch 폴더에 저장
stD = clDir(exp_idx, phantom_idx, acqType_idx, ang_deg, '..\Data_64ch\');
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

% Call function:: GenRFdata (128채널로 생성)
stRcvData = GenRFdata(stTrans, stTx, stScat, bPlot);

% 64채널 설정: 128채널에서 한 칸씩 건너뛰기 (1,3,5,...,127 또는 2,4,6,...,128)
nOrigChannel = 128;
nNewChannel = 64;
aChIdx = 1:2:nOrigChannel;  % 홀수 채널 선택 (1,3,5,...,127)
% aChIdx = 2:2:nOrigChannel;  % 짝수 채널 선택하려면 이 줄 사용 (2,4,6,...,128)

hFig = figure; hAx = axes;
% save mRcvData (64채널로 저장)
for txidx = 1:stRcvData.stTx.nNum
    mRcvData_128ch = stRcvData.vRcvData(:,:,txidx);
    
    % 64채널로 변환: 한 칸씩 건너뛰어 선택
    mRcvData = mRcvData_128ch(:, aChIdx);
    
    save([stD.sDataDir 'RcvData\' 'mRcvData_tx_' num2str(txidx) '.mat'],'mRcvData');
    imagesc(hAx, mRcvData); title(['mRcvData-tx-' num2str(txidx) ' (64ch), angle=' num2str(stTx.aPwAngle_deg(txidx)) ' deg']); 
    pause(0.1);
end

% Transducer 정보 업데이트 (64채널용)
stTrans_64ch = stTrans;
% 트랜스듀서 element 위치도 한 칸씩 건너뛰어 업데이트
if isfield(stTrans, 'aElemPos')
    stTrans_64ch.aElemPos = stTrans.aElemPos(aChIdx, :);
end
if isfield(stTrans, 'mElePos')
    stTrans_64ch.mElePos = stTrans.mElePos(aChIdx, :);
end
if isfield(stTrans, 'nElement')
    stTrans_64ch.nElement = nNewChannel;
end
% nNumEle, nNumEle_x 업데이트 (SetBFParam_simul에서 사용)
if isfield(stTrans, 'nNumEle')
    stTrans_64ch.nNumEle = nNewChannel;
end
if isfield(stTrans, 'nNumEle_x')
    stTrans_64ch.nNumEle_x = nNewChannel;
end
% pitch가 있다면 2배로 (한 칸 건너뛰므로)
if isfield(stTrans, 'nPitch')
    stTrans_64ch.nPitch = stTrans.nPitch * 2;
end
if isfield(stTrans, 'nPitch_x')
    stTrans_64ch.nPitch_x = stTrans.nPitch_x * 2;
end

% save stRFInfo (64채널 정보로 저장)
stRFInfo.sVersion       = stRcvData.sVersion;
stRFInfo.stTrans        = stTrans_64ch;  % 64채널 트랜스듀서 정보
stRFInfo.nRcvChannel    = nNewChannel;   % 64채널로 변경
stRFInfo.nSoundSpeed    = stRcvData.nSoundSpeed;
stRFInfo.nSample        = size(mRcvData,1);
stRFInfo.nFs            = stRcvData.nFs;
stRFInfo.nFc            = stRcvData.nFc;

% stTx 정보 복사 및 64채널로 변환
stTx_64ch = stRcvData.stTx;
% mAptr (Tx aperture) 64채널로 변환: [nNum x 128] -> [nNum x 64]
if isfield(stRcvData.stTx, 'mAptr')
    stTx_64ch.mAptr = stRcvData.stTx.mAptr(:, aChIdx);
end
% mDelay (Tx delay) 64채널로 변환: [nNum x 128] -> [nNum x 64]
if isfield(stRcvData.stTx, 'mDelay')
    stTx_64ch.mDelay = stRcvData.stTx.mDelay(:, aChIdx);
end
stTx_64ch.aRxAptrPos = zeros(1, stRcvData.stTx.nNum);  % Rx aperture position (64ch에서 필요)

stRFInfo.stTx           = stTx_64ch;
stRFInfo.stScat         = stRcvData.stScat;
stRFInfo.sStartAcqTime  = stRcvData.sStartAcqTime;
stRFInfo.aChIdx         = aChIdx;        % 사용된 채널 인덱스 저장 (참고용)
save([stD.sDataDir 'RcvData\' 'stRFInfo.mat'],'stRFInfo');

disp(['64-channel data saved to: ' stD.sDataDir]);
disp(['Selected channels: ' num2str(aChIdx(1)) ':2:' num2str(aChIdx(end))]);
