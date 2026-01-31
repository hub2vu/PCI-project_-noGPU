% Sua Bae
%  2019-05-22
%     BFpm:  Ver 3.0
%  2025-01-31
%     Modified for 64 channel data
% 
clc;clear; close all;
addpath('src');
addpath('src\beamformer');

%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nMultiBeam = 1;
        
% transducer
load('stTrans_L7-4.mat');

exp_idx = 1;        % 1:simulation, 2:phantom, 3:invivo
phantom_idx = 3;
acqType_idx = 2;    % 1:CF, 2:PW
ang_deg = 4;        % 0 or 10 (deg of max PW ang)
             
% directory - 64ch 데이터 경로
stD = clDir(exp_idx, phantom_idx, acqType_idx, ang_deg, '..\Data_64ch');
CreateFolders(stD);

cBF = {'RFBF','IQBF'};

bf_idx = 2; %  1:RFBF, 2:IQBF
nQiq = 2;
nQph = 0;

sRcvDir = [stD.sDataDir 'RcvData\'];
sBfDir  = [stD.sDataDir 'BfData\' cBF{bf_idx} '\'];

sBfType = cBF{bf_idx};

if ~exist(sBfDir,'dir');mkdir(sBfDir);end

% load RFInfo (64채널 정보)
load([sRcvDir 'stRFInfo.mat']);

% 64채널 확인 메시지
disp(['Processing ' num2str(stRFInfo.nRcvChannel) ' channel data...']);

% call function
if strcmp(sBfType,'RFBF')
    disp('RF beamforming ....');
    stBFpm = SetBFParam_simul(sBfType, stRFInfo, nMultiBeam, nQiq, nQph);
    RFBeamformer(sRcvDir, sBfDir, stRFInfo, stBFpm);
elseif strcmp(sBfType,'IQBF')
    disp('IQ beamforming ....');
    stBFpm = SetBFParam_simul(sBfType, stRFInfo, nMultiBeam, nQiq, nQph);
    IQBeamformer(sRcvDir, sBfDir, stRFInfo, stBFpm);
end

% Load
load([stD.sDataDir 'BfData\' stBFpm.sBfType '\' stBFpm.sBfOutName]);
mBfData = stBfData.mBfData;
aX = stBfData.stBFpm.stG.aX;
aZ = stBfData.stBFpm.stG.aZ;   


% normalize
mBfData_norm = abs(mBfData)/max(mBfData(:));

% plot image
figure('color','w');
imagesc(aX*1e3, aZ*1e3, 20*log10(mBfData_norm));
clim([-60 0]); colormap(gray); axis equal; axis tight;
set(gca,'FontName','Times New Roman', 'FontSize',12);
xlabel('lateral position (mm)')
ylabel('axial position (depth) (mm)')
title(['64-channel Beamformed Image (' sBfType ')']);
