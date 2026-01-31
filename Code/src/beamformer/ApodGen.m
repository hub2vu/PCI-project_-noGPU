function aApod = ApodGen(aDistance,aApertureSize,sWindowType)

    switch(sWindowType)
        case 'none' 
            aApod = ones(size(aDistance)); 
        case 'boxcar' 
            aApod = double(aDistance<=aApertureSize/2); 
        case 'hanning'
            aApod = double(aDistance<=aApertureSize/2).*(0.5 + 0.5*cos(2*pi*aDistance./aApertureSize)); 
        case 'hamming'
            aApod = double(aDistance<=aApertureSize/2).*(0.53836 + 0.46164*cos(2*pi*aDistance./aApertureSize)); 
        case 'tukey25'
            roll=0.25;
            aApod =(aDistance<(aApertureSize/2*(1-roll))) + (aDistance>(aApertureSize/2*(1-roll))).*(aDistance<(aApertureSize/2)).*0.5.*(1+cos(2*pi/roll*(aDistance./aApertureSize-roll/2-1/2)));                               
        case 'tukey50'
            roll=0.5;
            aApod=(aDistance<(aApertureSize/2*(1-roll))) + (aDistance>(aApertureSize/2*(1-roll))).*(aDistance<(aApertureSize/2)).*0.5.*(1+cos(2*pi/roll*(aDistance./aApertureSize-roll/2-1/2)));                               
        case 'tukey75'
            roll=0.75;
            aApod=(aDistance<(aApertureSize/2*(1-roll))) + (aDistance>(aApertureSize/2*(1-roll))).*(aDistance<(aApertureSize/2)).*0.5.*(1+cos(2*pi/roll*(aDistance./aApertureSize-roll/2-1/2)));                               
        case 'plancktaper'
            e = 0.004./aApertureSize;
%             e = 0.003./aApertureSize;
%             e = 0.002./aApertureSize;
            z = 2*e.*( 1./(1-2*aDistance./aApertureSize) + 1./(1-2*e-2*aDistance./(aApertureSize)) );
            aApod = (aDistance<=((0.5-e).*aApertureSize)) + (aDistance>((0.5-e).*aApertureSize)).*(aDistance<(aApertureSize/2)).*(1./(exp(z)+1));
        otherwise
            error('Unknown window type. Known types are: boxcar, hamming, hanning, tukey25, tukey50, tukey75.');
    end
    
end
