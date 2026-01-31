"""
UNet_DataLoader.py

MATLAB에서 생성된 128ch/64ch 데이터셋을 PyTorch U-Net 학습에 사용하기 위한 데이터 로더

사용법:
    from UNet_DataLoader import UltrasoundDataset
    
    dataset = UltrasoundDataset(data_dir='../dataset/')
    train_loader = DataLoader(dataset, batch_size=8, shuffle=True)

2025-01-31
"""

import os
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from scipy.io import loadmat
import cv2


class UltrasoundDataset(Dataset):
    """
    초음파 64ch -> 128ch 변환을 위한 데이터셋
    
    Args:
        data_dir: 데이터셋 루트 디렉토리 (bf_128ch/, bf_64ch/, metadata/ 포함)
        transform: 데이터 증강 변환 (optional)
        target_size: 출력 이미지 크기 (H, W) - None이면 원본 크기 유지
        normalize: 정규화 방식 ('minmax', 'zscore', None)
    """
    
    def __init__(self, data_dir, transform=None, target_size=None, normalize='minmax'):
        self.data_dir = data_dir
        self.transform = transform
        self.target_size = target_size
        self.normalize = normalize
        
        # 데이터 파일 목록 가져오기
        self.bf_128ch_dir = os.path.join(data_dir, 'bf_128ch')
        self.bf_64ch_dir = os.path.join(data_dir, 'bf_64ch')
        self.metadata_dir = os.path.join(data_dir, 'metadata')
        
        # 파일 ID 목록
        self.file_ids = []
        for f in os.listdir(self.bf_128ch_dir):
            if f.endswith('.mat'):
                # bf_000001.mat -> 000001
                file_id = f.replace('bf_', '').replace('.mat', '')
                self.file_ids.append(file_id)
        
        self.file_ids.sort()
        print(f"Found {len(self.file_ids)} samples in {data_dir}")
        
    def __len__(self):
        return len(self.file_ids)
    
    def __getitem__(self, idx):
        file_id = self.file_ids[idx]
        
        # 128ch (Ground Truth) 로드
        mat_128ch = loadmat(os.path.join(self.bf_128ch_dir, f'bf_{file_id}.mat'))
        img_128ch = mat_128ch['mBfData_128ch_norm'].astype(np.float32)
        
        # 64ch (Input) 로드
        mat_64ch = loadmat(os.path.join(self.bf_64ch_dir, f'bf_{file_id}.mat'))
        img_64ch = mat_64ch['mBfData_64ch_norm'].astype(np.float32)
        
        # 크기 조정 (64ch를 128ch 크기에 맞춤)
        if img_64ch.shape != img_128ch.shape:
            img_64ch = cv2.resize(img_64ch, (img_128ch.shape[1], img_128ch.shape[0]), 
                                   interpolation=cv2.INTER_LINEAR)
        
        # 타겟 크기로 리사이즈
        if self.target_size is not None:
            img_128ch = cv2.resize(img_128ch, (self.target_size[1], self.target_size[0]),
                                    interpolation=cv2.INTER_LINEAR)
            img_64ch = cv2.resize(img_64ch, (self.target_size[1], self.target_size[0]),
                                   interpolation=cv2.INTER_LINEAR)
        
        # 정규화
        if self.normalize == 'minmax':
            # 이미 0~1 범위이지만 확인
            img_128ch = np.clip(img_128ch, 0, 1)
            img_64ch = np.clip(img_64ch, 0, 1)
        elif self.normalize == 'zscore':
            img_128ch = (img_128ch - img_128ch.mean()) / (img_128ch.std() + 1e-8)
            img_64ch = (img_64ch - img_64ch.mean()) / (img_64ch.std() + 1e-8)
        
        # 채널 차원 추가 [H, W] -> [1, H, W]
        img_128ch = np.expand_dims(img_128ch, axis=0)
        img_64ch = np.expand_dims(img_64ch, axis=0)
        
        # 텐서 변환
        input_tensor = torch.from_numpy(img_64ch)
        target_tensor = torch.from_numpy(img_128ch)
        
        # 데이터 증강 적용
        if self.transform is not None:
            # 동일한 변환을 input과 target에 적용
            seed = np.random.randint(2147483647)
            torch.manual_seed(seed)
            input_tensor = self.transform(input_tensor)
            torch.manual_seed(seed)
            target_tensor = self.transform(target_tensor)
        
        return {
            'input': input_tensor,      # 64ch (저해상도)
            'target': target_tensor,    # 128ch (고해상도)
            'file_id': file_id
        }


class UltrasoundDatasetWithMetadata(UltrasoundDataset):
    """
    메타데이터도 함께 반환하는 데이터셋
    """
    
    def __getitem__(self, idx):
        sample = super().__getitem__(idx)
        file_id = sample['file_id']
        
        # 메타데이터 로드
        mat_meta = loadmat(os.path.join(self.metadata_dir, f'meta_{file_id}.mat'))
        metadata = mat_meta['stMetadata']
        
        sample['num_scatterers'] = int(metadata['nNumScatterers'][0,0][0,0])
        sample['scatterer_positions'] = metadata['mScatPos_mm'][0,0]
        
        return sample


def create_dataloaders(data_dir, batch_size=8, train_ratio=0.8, target_size=(512, 256), 
                       num_workers=4):
    """
    학습/검증 데이터 로더 생성
    
    Args:
        data_dir: 데이터셋 디렉토리
        batch_size: 배치 크기
        train_ratio: 학습 데이터 비율
        target_size: (H, W) 이미지 크기
        num_workers: 데이터 로딩 워커 수
        
    Returns:
        train_loader, val_loader
    """
    
    # 전체 데이터셋 생성
    full_dataset = UltrasoundDataset(data_dir, target_size=target_size)
    
    # 학습/검증 분할
    n_total = len(full_dataset)
    n_train = int(n_total * train_ratio)
    n_val = n_total - n_train
    
    train_dataset, val_dataset = torch.utils.data.random_split(
        full_dataset, [n_train, n_val],
        generator=torch.Generator().manual_seed(42)
    )
    
    # 데이터 로더 생성
    train_loader = DataLoader(
        train_dataset, 
        batch_size=batch_size, 
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True
    )
    
    val_loader = DataLoader(
        val_dataset, 
        batch_size=batch_size, 
        shuffle=False,
        num_workers=num_workers,
        pin_memory=True
    )
    
    print(f"Train samples: {n_train}, Validation samples: {n_val}")
    
    return train_loader, val_loader


# 테스트
if __name__ == '__main__':
    # 데이터셋 테스트
    data_dir = '../dataset/'
    
    if os.path.exists(data_dir):
        dataset = UltrasoundDataset(data_dir, target_size=(512, 256))
        
        # 첫 번째 샘플 확인
        sample = dataset[0]
        print(f"Input shape: {sample['input'].shape}")
        print(f"Target shape: {sample['target'].shape}")
        print(f"Input range: [{sample['input'].min():.4f}, {sample['input'].max():.4f}]")
        print(f"Target range: [{sample['target'].min():.4f}, {sample['target'].max():.4f}]")
        
        # 데이터 로더 테스트
        train_loader, val_loader = create_dataloaders(data_dir, batch_size=4)
        
        for batch in train_loader:
            print(f"Batch input shape: {batch['input'].shape}")
            print(f"Batch target shape: {batch['target'].shape}")
            break
    else:
        print(f"Data directory not found: {data_dir}")
        print("Run Generate_Training_Dataset.m first to create the dataset.")
