% 1. 데이터 불러오기
% (nrrdread 함수가 같은 폴더에 있어야 합니다)
filename = 'annotation_25.nrrd';
data = nrrdread(filename);

% 2. 데이터 크기 줄이기 (다운샘플링)
% 원본(528x320x456)은 너무 커서 3D로 돌릴 때 렉이 걸립니다.
% 4분의 1 크기로 줄여서 가볍게 확인합니다.
scale_factor = 0.25; 
vol_small = imresize3(double(data), scale_factor, 'nearest');

% 배경(0)과 뇌 조직(0보다 큰 값) 구분
brain_mask = vol_small > 0;

%% 방법 1: Volshow 사용 (가장 추천)
% MATLAB 2017a 이상 Image Processing Toolbox 필요
try
    figure;
    % 뇌 영역을 입체적으로 렌더링
    volshow(brain_mask, 'Renderer', 'VolumeRendering', 'BackgroundColor', [0 0 0]);
    title('Allen Mouse Brain Atlas (3D Volume)');
catch
    disp('volshow 함수가 없습니다. 방법 2로 진행합니다.');
end

%% 방법 2: Isosurface 사용 (모든 버전 호환)
% 뇌의 겉 표면(껍데기)을 그려줍니다.
figure;
p = patch(isosurface(brain_mask, 0.5));
p.FaceColor = [0.8 0.4 0.4]; % 붉은색 뇌
p.EdgeColor = 'none';
p.FaceAlpha = 0.6; % 약간 투명하게

% 조명 및 카메라 설정
daspect([1 1 1]);
view(3); 
axis tight;
camlight left; 
lighting gouraud;
grid on;
title('Mouse Brain Surface');
xlabel('X (Lateral)'); ylabel('Y (Vertical)'); zlabel('Z (Anterior)');

% 마우스로 회전시키며 확인해보세요 (Rotate 3D 아이콘 클릭)
rotate3d on;