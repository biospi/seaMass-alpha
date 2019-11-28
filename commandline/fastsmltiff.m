function [x,y,img,ast] = fastsmltiff(filename)
    % fastsmltiff - Load sml File generated from a TIFF fike.
    %
    % Syntax: [x,y,img,ast] = fastsmltiff(filename)
    %
    % Output image index for x and y, raw image - img - and 
    % Anscombe Transform - ast - log(2*sqrt(sc'+3/8)))
    % Plot using:
    % 	imagesc(img);
    % 	imagesc(ast); set(gca,'YDir','normal')

    
    specIdx='/mzScale=0/spectrumIndex';
    binLoc='/mzScale=0/binLocations';
    binInt='/mzScale=0/binIntensities';
    
    scanIdx=h5read(filename,specIdx);
    xdata=h5read(filename,binLoc);
    ydata=h5read(filename,binInt);
    
    minX = min(xdata)+1;
    maxX = max(xdata)+1;
    
    clear x;
    
    m = length(scanIdx)-1;
    n = maxX;
    
    img = zeros(m,n);
    
    x=[1:1:n];
    y=[1:1:m];
    
    for i = y;
        dataPos = scanIdx(i);
        dataLen = scanIdx(i+1);
        %idx = h5read(filename,binLoc,double(dataPos+1),double(dataLen));
        %inten = h5read(filename,binInt,double(dataPos+1),double(dataLen));
        idx = xdata(dataPos+1:dataLen);
        inten = ydata(dataPos+1:dataLen);
        disp(['Reading scan: ',num2str(i)]);
        if length(idx) < maxX
            disp(['Rebinning scan: ',num2str(i)]);
            intIdx = 1;
            for j = 1:length(idx)
                img(i,idx(j)+1) = inten(intIdx);
                intIdx = intIdx + 1;
            end
        else
           img(i,:)=inten; 
        end
    end
    ast = log(2.*sqrt(img+3/8));
end
