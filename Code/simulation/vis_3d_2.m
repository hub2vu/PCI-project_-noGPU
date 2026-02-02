% 1) 데이터 로드
data = nrrdread('annotation_25.nrrd');
sz = size(data);

midX = round(sz(1)/2);
midZ = round(sz(3)/2);

% 2) 세로 단면(예: sagittal)
sl_v = data(:,:,midZ);

% 3) 가로 단면(예: axial)  -> X축에서 자른 단면
sl_h = squeeze(data(midX,:,:))';   % ' 를 붙여서 보기 좋게 회전(필요 없으면 빼도 됨)

% 4) 경계 계산 함수(간단 inline)
getBoundary = @(sl) ...
    ( ...
      ( [sl(:,1:end-1) ~= sl(:,2:end), false(size(sl,1),1)] & (sl>0) ) | ...
      ( [sl(1:end-1,:) ~= sl(2:end,:); false(1,size(sl,2))] & (sl>0) ) ...
    );

Bv = getBoundary(sl_v);
Bh = getBoundary(sl_h);

% 5) 시각화
figure;

% --- 세로
subplot(1,2,1);
mask_v = sl_v > 0;
imshow(mask_v); hold on;
[y,x] = find(Bv);
plot(x,y,'r.','MarkerSize',1);
title('Vertical slice + boundaries');

% --- 가로
subplot(1,2,2);
mask_h = sl_h > 0;
imshow(mask_h); hold on;
[y,x] = find(Bh);
plot(x,y,'r.','MarkerSize',1);
title('Horizontal slice + boundaries');
