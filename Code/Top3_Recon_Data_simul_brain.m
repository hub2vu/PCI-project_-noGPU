% Top3_Recon_Data_simul_brain.m
% Load RF -> SetBFParam_simul -> ROI 제한 -> IQ Beamforming -> B-mode plot

clc; clear; close all;

addpath('src');
addpath('src\beamformer');

%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nMultiBeam = 1;

% transducer
load('stTrans_L7-4.mat');

exp_idx = 1;        % 1:simulation
phantom_idx = 3;    % Top1과 동일하게 유지
acqType_idx = 2;    % 1:CF, 2:PW
ang_deg = 4;        % Top1과 동일 필수

% directory (Top1과 동일 루트)
stD = clDir(exp_idx, phantom_idx, acqType_idx, ang_deg, '..\Data_brain');
CreateFolders(stD);

cBF = {'RFBF','IQBF'};
bf_idx = 2; % 1:RFBF, 2:IQBF
nQiq = 2;
nQph = 0;

sRcvDir = [stD.sDataDir 'RcvData\'];
sBfDir  = [stD.sDataDir 'BfData\' cBF{bf_idx} '\'];
sBfType = cBF{bf_idx};

if ~exist(sBfDir,'dir'); mkdir(sBfDir); end

% load RFInfo
load([sRcvDir 'stRFInfo.mat']);

% ===== Beamforming Params =====
stBFpm = SetBFParam_simul(sBfType, stRFInfo, nMultiBeam, nQiq, nQph);

% ===== ROI 제한(중요: 속도/메모리) =====
X_MAX_MM = 7.0;     % lateral +/- 7mm
Z_MAX_MM = 12.0;    % depth 0~12mm

x_max = X_MAX_MM * 1e-3;
z_max = Z_MAX_MM * 1e-3;

dx = stBFpm.stG.dx;
dz = stBFpm.stG.dz;

stBFpm.stG.aX = (-x_max:dx:x_max);
stBFpm.stG.nXdim = numel(stBFpm.stG.aX);

stBFpm.stG.aZ = (0:dz:z_max);
stBFpm.stG.nZdim = numel(stBFpm.stG.aZ);

disp("IQ beamforming ....");
IQBeamformer(sRcvDir, sBfDir, stRFInfo, stBFpm);

% Load BF output
load([stD.sDataDir 'BfData\' stBFpm.sBfType '\' stBFpm.sBfOutName]);
mBfData = stBfData.mBfData;
aX = stBfData.stBFpm.stG.aX;
aZ = stBfData.stBFpm.stG.aZ;

% normalize + B-mode
mBfData_norm = abs(mBfData) / (max(abs(mBfData(:))) + eps);

figure('color','w');
imagesc(aX*1e3, aZ*1e3, 20*log10(mBfData_norm + 1e-12));
clim([-60 0]); colormap(gray); axis equal; axis tight;
set(gca,'FontName','Times New Roman', 'FontSize',12);
xlabel('lateral position (mm)');
ylabel('axial position (depth) (mm)');
title('Brain phantom B-mode');
