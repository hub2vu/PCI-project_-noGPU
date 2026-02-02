clear; close all; clc;

disp('=== START verify_brain_phantom ===');
disp(['pwd = ', pwd]);

% NOTE: 상위폴더 전체 genpath는 느릴 수 있어 필요 최소만 추가 권장
addpath(pwd);

annPath   = fullfile(pwd, 'annotation_25.nrrd');
cachePath = fullfile(pwd, 'brain_phantom_fieldii.mat');

disp(['annPath   = ', annPath]);
disp(['cachePath = ', cachePath]);

if ~isfile(annPath)
    error(['annotation_25.nrrd not found: ', annPath, ' (현재 폴더가 simulation인지 확인)']);
end

if isfile(cachePath)
    disp('Cache exists. Checking variables inside...');
    whos('-file', cachePath);

    disp('Loading stScat/meta...');
    S = load(cachePath);
    if ~isfield(S,'stScat') || ~isfield(S,'meta')
        error('Cache file exists but does not contain stScat/meta. Delete cache and rebuild.');
    end
    stScat = S.stScat;
    meta   = S.meta;
else
    disp('Cache NOT found. Building phantom...');

    % build 함수가 있는지 확인
    w = which('build_brain_phantom_fieldii');
    if isempty(w)
        error('build_brain_phantom_fieldii.m 을 찾을 수 없습니다. simulation 폴더에 파일이 있는지 확인하세요.');
    else
        disp(['found build_brain_phantom_fieldii: ', w]);
    end

    opts.Nscat = 200000;
    opts.vox_um = 25;
    opts.slab_mm = 0.4;
    opts.z0_mm = 0.6;
    opts.amp_scale = 1e26;
    opts.amp_sigma = 0.6;
    opts.seed = 0;

    [stScat, meta] = build_brain_phantom_fieldii(annPath, cachePath, opts);
    disp(['Saved phantom cache: ', cachePath]);
end

mask = meta.mask;
disp(['mask size = ', mat2str(size(mask))]);
disp(['mask nnz  = ', num2str(nnz(mask))]);

if nnz(mask) == 0
    error('mask가 비었습니다. slab_mm을 키우거나 nrrd 읽기/경로 문제를 확인하세요.');
end

% 1) mask 단면 확인
Ny = size(mask,2);
y0 = round(Ny/2);
mask_xz = squeeze(mask(:, y0, :)); % [Nx x Nz]

figure('Color','w','Name','mask_xz');
set(gcf,'Position',[100 100 900 700]);
imagesc(mask_xz'); axis image; colormap(gray);
xlabel('x (voxel)'); ylabel('z (voxel)'); title('Brain mask slice (XZ)');
drawnow;

% 2) 산란체 3D 확인
xyz = stScat.mScatXYZPos;
disp(['xyz size = ', mat2str(size(xyz))]);

x = xyz(:,1)*1e3;
y = xyz(:,2)*1e3;
z = xyz(:,3)*1e3;

figure('Color','w','Name','scatter3');
set(gcf,'Position',[1050 100 900 700]);
scatter3(x, y, z, 1, z, 'filled');
axis equal; grid on; view(3);
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
title('Field II scatterer cloud (brain-shaped)');
drawnow;

fprintf('scatterers=%d, z(mm)=%.3f~%.3f\n', size(xyz,1), min(z), max(z));
disp('=== DONE ===');
