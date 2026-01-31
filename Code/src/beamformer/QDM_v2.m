% Sua Bae
% 2019-04-05
%    short-time FT function is now nested
%    In case of fixed demodulation, fftshift and fft are uesd
% 2019-02-07
%    Demodulation with a fixed freq OR Dynamic Demodulation with a linear decay
%
function [mOutData] = QDM_v2(mData, nFs, nSoundSpeed, stDemod, nDeci)

    bPlot = 0;

    if stDemod.bDynamic
        nDemodFreq_slop         = stDemod.nDemodFreq_slop; % [Hz/m]
        nDemodFreq_intercept    = stDemod.nDemodFreq_intercept; % [Hz]
        nDemodLPFFreq           = stDemod.nDemodLPFFreq;
        nLPFTap            = stDemod.nDemodLPFTap;
        
        %% Demodulation    
%         nPhase = 0;
%         for n = 0:1:size(mData,1)-1; 
%             nDepth = n/nFs*nSoundSpeed/2; % [m]
%             nDemodFreq = nDemodFreq_slop*nDepth + nDemodFreq_intercept; % [Hz]
%             nPhase = nPhase + 2*pi*nDemodFreq/nFs; % [rad]
%             cos_n(n+1) = cos(nPhase);
%             sin_n(n+1) = sin(nPhase);
%         end
        aK = (1:size(mData,1));
        aPhase = pi*nDemodFreq_slop*nSoundSpeed/nFs/nFs*aK.*(aK-1)/2 + aK*2*pi*nDemodFreq_intercept/nFs;
        cos_n = cos(aPhase);
        sin_n = sin(aPhase);
        
        mInph_temp       = bsxfun(@times,mData, cos_n');
        mQuad_temp    = bsxfun(@times,mData, -sin_n');
        
        if(bPlot)
            figure('position',[200 200 800 500]);
            subplot(2,2,1);ShortTimeFT(mData(:,32), 512, 256, nFs);title('Original data');
                           hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nDemodFreq_intercept,aK,'r');
                           hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nDemodFreq_intercept - (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
                           hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nDemodFreq_intercept + (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
            subplot(2,2,2);ShortTimeFT(mInph_temp(:,32), 512, 256, nFs);title('after dynamic demod');
                           hold on; plot(- (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
                           hold on; plot(+ (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
        end
        
        if 0
            % plot 3 freq bands for Freq Compounding 
            figure('position',[200 200 1200 300]);
            nintrcpt = nDemodFreq_intercept;
            subplot(1,3,1);ShortTimeFT(mData(:,32), 512, 256, nFs);title('Original data');
            hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nintrcpt,aK,'r');
            hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nintrcpt - (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
            hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nintrcpt + (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
            xlim([-5 5]*1e6);
            nintrcpt = nintrcpt - 0.5e6;
            subplot(1,3,2);ShortTimeFT(mData(:,32), 512, 256, nFs);title('Original data');
            hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nintrcpt,aK,'r');
            hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nintrcpt - (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
            hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nintrcpt + (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
            xlim([-5 5]*1e6);
            nintrcpt = nintrcpt +1e6;
            subplot(1,3,3);ShortTimeFT(mData(:,32), 512, 256, nFs);title('Original data');
            hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nintrcpt,aK,'r');
            hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nintrcpt - (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
            hold on; plot(nDemodFreq_slop*aK/nFs*nSoundSpeed/2+nintrcpt + (nDemodLPFFreq + nDemodFreq_slop/2*aK/nFs*nSoundSpeed/2),aK,'k');
            xlim([-5 5]*1e6);
        end
        %% dynamic LPF & Decimation 
        % filter len
        nLPFTap = floor(nLPFTap/2)*2+1; % odd integer!
        nTap_half = floor(nLPFTap/2);
        % input (zero-padding)
        mInph_zp = [zeros(nTap_half,size(mData,2));mInph_temp;zeros(nTap_half,size(mData,2))];
        mQuad_zp = [zeros(nTap_half,size(mData,2));mQuad_temp;zeros(nTap_half,size(mData,2))];
        % output
        nLen_deci = floor(size(mData,1)/nDeci);
        mInph   = zeros(nLen_deci,size(mData,2));
        mQuad   = zeros(nLen_deci,size(mData,2));
        % filtering and decimation
        % Cutoff freq of LPF is decaying by a factor of slop/2
        for oidx = 1:nLen_deci % decimated idx
            nLPFFreq = nDemodLPFFreq + nDemodFreq_slop/1.8*(oidx-1)*nDeci/nFs*nSoundSpeed/2; 
            aLPFCoefs = fir1(nLPFTap-1, nLPFFreq/nFs*2, 'low')';
            sidx = (oidx-1)*nDeci+1;
            mInph(oidx,:) = sum(bsxfun(@times,mInph_zp(sidx:sidx+nLPFTap-1,:),aLPFCoefs),1);
            mQuad(oidx,:) = sum(bsxfun(@times,mQuad_zp(sidx:sidx+nLPFTap-1,:),aLPFCoefs),1);
        end
%         aLPFCoefs = fir1(nLPFTap-1, nDemodLPFFreq/nFs*2, 'low');
%         mInph = convn(mInph_temp, aLPFCoefs', 'same');
%         mQuad = convn(mQuad_temp, aLPFCoefs', 'same');
        %% Complex
        mOutData = mInph + 1i*mQuad;            
        
        if(bPlot)
            subplot(2,2,3);ShortTimeFT(mInph(:,32), round(512/nDeci), round(256/nDeci), nFs/nDeci);  title('after dynamic LPF and deci');
            subplot(2,2,4);ShortTimeFT(mOutData(:,32), round(512/nDeci), round(256/nDeci), nFs/nDeci); title('I +jQ');
        end
    else
        nDemodFreq = stDemod.nDemodFreq; % [Hz]
        nDemodLPFFreq = stDemod.nDemodLPFFreq;
        nLPFTap = stDemod.nDemodLPFTap;
        
        %% Demodulation    
        n =0:size(mData,1)-1; 
        cos_n = cos(2*pi*nDemodFreq/nFs*n);
        sin_n = sin(2*pi*nDemodFreq/nFs*n);
        
        mInph_temp       = bsxfun(@times,mData, cos_n');
        mQuad_temp    = bsxfun(@times,mData, -sin_n');
        
        %% LPF    
        aLPFCoefs = fir1(nLPFTap-1, nDemodLPFFreq/nFs*2, 'low');
        mInph = convn(mInph_temp, aLPFCoefs', 'same');
        mQuad = convn(mQuad_temp, aLPFCoefs', 'same');

        %% Decimation 
        mOutData = mInph(1:nDeci:end, :) + 1i*mQuad(1:nDeci:end, :); 
        
        if bPlot
            hFig = figure(99);
            set(hFig, 'Position', [300 300 1000 400]);
            subplot(2,3,1);
                plot(linspace(-nFs/2,nFs/2,4096)*1e-6,db(abs(fftshift(fft(mData(:,64),4096))))); nPlotMax = max(db(abs(fftshift(fft(mData(:,64),4096)))));
                xlim([-nFs/2,nFs/2]);xlim([-nFs/2,nFs/2]*1e-6); ylim([nPlotMax-120 nPlotMax]);grid on; grid minor;title('original data');
            subplot(2,3,2);
                plot(linspace(-nFs/2,nFs/2,4096)*1e-6,db(abs(fftshift(fft((mInph_temp(:,64)+1j*mQuad_temp(:,64)),4096)))));xlim([-nFs/2,nFs/2]*1e-6); ylim([nPlotMax-120 nPlotMax]);;grid on; grid minor;title('demodulated');
            subplot(2,3,4);
                plot(linspace(-nFs/2,nFs/2,4096)*1e-6,db(abs(fftshift(fft(aLPFCoefs,4096)))));xlim([-nFs/2,nFs/2]);xlim([-nFs/2,nFs/2]*1e-6);grid on; grid minor;title('filter');
            subplot(2,3,5);
                plot(linspace(-nFs/2,nFs/2,4096)*1e-6,db(abs(fftshift(fft((mInph(:,64)+1j*mQuad(:,64)),4096)))));xlim([-nFs/2,nFs/2]*1e-6); ylim([nPlotMax-120 nPlotMax]);grid on; grid minor;title('filtered');
            subplot(2,3,6);
                plot(linspace(-nFs/2/nDeci,nFs/2/nDeci,4096)*1e-6,db(abs(fftshift(fft(mOutData(:,64),4096)))));xlim([-nFs/2/nDeci,nFs/2/nDeci]);xlim([-nFs/2/nDeci,nFs/2/nDeci]*1e-6); ylim([nPlotMax-120 nPlotMax]);grid on; grid minor;title('filtered+decimated');
            % figure;
            % subplot(2,2,1);ShortTimeFT(mData(:,32), 512, 256, nFs);title('Original data');
            % subplot(2,2,2);ShortTimeFT(mInph_temp(:,32), 512, 256, nFs);title('after fixed demod');
            % subplot(2,2,3);ShortTimeFT(mInph(:,32), round(512), round(256), nFs);  title('after fixed LPF');
            % subplot(2,2,4);ShortTimeFT(mOutData(:,32), round(512/nDeci), round(256/nDeci), nFs/nDeci); title('I +jQ and deci')
        end
    end


%     function ShortTimeFT(aData, nFftLen, nShiftLen, nFs)
%     %     figure;
%         nDataLen = numel(aData);
%         nOutLen = ceil(nDataLen/nShiftLen);
%         for tidx = 1:nOutLen
%             sidx = (tidx-1)*nShiftLen+1;
%             eidx = min(sidx+nFftLen-1,nDataLen);
%             mSTFT(tidx,:) = fftshift(abs(fft(aData(sidx:eidx),2048)));
%             aTAx(tidx) = ( sidx + eidx )/2;
%         end
%     %     subplot(1,2,1);
%     %         imagesc(linspace(-0.5,0.5,nFftLen)*nFs,aTAx,db(mSTFT));
%     %         ylabel('sample'); xlabel('frequency'); colorbar;
%     %     subplot(1,2,2);
%             mSTFT_norm = bsxfun(@times,mSTFT,1./max(mSTFT,[],2));
%             imagesc(linspace(-0.5,0.5,nFftLen)*nFs,aTAx,db(mSTFT_norm));
%             ylabel('sample'); xlabel('frequency'); colorbar; caxis([-30 0]);
% 
%     end
    
end