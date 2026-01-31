% Sua Bae
%
% 2019-05-22
%   update: bypass
% 2019-04-05
%   Interpolation
%
function [mOutData] = Interp_v1(mData, nIntpFactor, nFiltTap)

    if ~exist('nFiltTap','var')
        nFiltTap = 257;
    end

    if nIntpFactor == 1 % bypass
        mOutData = mData;
        
    else
        
        mData_padded = zeros(size(mData,1)*nIntpFactor, size(mData,2));
        mData_padded(1:nIntpFactor:end,:) = mData;

        aLPFCoefs = fir1(nFiltTap-1, 1/nIntpFactor, 'low');
        mOutData = convn(mData_padded, aLPFCoefs', 'same') *nIntpFactor; % multiplying 'nIntpFactor' for maintaining amplitude

        bPlot = 0;
        if bPlot
            figure;
                plot(real(mData_padded(:,round(size(mData,2)/2))));
                hold on; plot(real(mOutData(:,round(size(mData,2)/2)))); grid on; grid minor;
            figure;
                freqz(mData_padded(:,round(size(mData,2)/2)),1,4096);
            figure;
                freqz(mOutData(:,round(size(mData,2)/2)),1,4096);
        end
    end
    
end