% Sua Bae
%  2019-05-22
% 
clc;clear; close all;
addpath('src');


% directory
exp_idx = 1;        % 1:simulation, 2:phantom, 3:invivo
phantom_idx = 3;
acqType_idx = 2;    % 1:CF, 2:PW
ang_deg = 4;        % 0 or 10 (deg of max PW ang)

sBFType = 'IQBF';
nBW = 80;
aQiq = 2;

stD     = clDir(exp_idx, phantom_idx, acqType_idx, ang_deg, '..\Data');   
sBfDir  = [stD.sDataDir 'BfData\' sBFType '\'];
sBfOutName = 'stBfData.mat';


% load data
load([sBfDir sBfOutName]);
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
