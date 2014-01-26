function [ synthFade,uvFade,outFadePos ] = bvFade( inFadePos,wLen,hopSize,bvVUVbool,lvVUVbool,fadeLength,expon )
%MATTFADE  Calculates crossfade arrays for transitions between voiced and
%unvoiced sections. 
% inFadePos: number between 0 and 1 specifying fade position. 1 is all
% unvoiced unsynthesized waveform. 0 is all synthesized.
% wLen: length of synthesis window
% hopSize: the hopsize specifies how long the output vectors will be
% bvVUVbool/lvVUVbool: booleans specifying whether current BV or LV windows
% are voiced.
% fadelength: the length over which the crossfade occurs. Specified in
% amount of windows (ie 1-3). Must be a whole number
% expon: exponent controling the curve of the crossfade. 1 is a linear
% crossfade. 2 produces good results.

synthFade = zeros(hopSize,1);
uvFade = zeros(hopSize,1);
fadeOutTable = (1/(fadeLength*wLen)):(1/(fadeLength*wLen)):1;
fadeOutTable = 1-power(fadeOutTable,expon);
fadeInTable = (1-fadeOutTable);


if (bvVUVbool==0 || lvVUVbool==0)% if either unvoiced
           endPos = inFadePos-1+hopSize;
    if(endPos < wLen*fadeLength-1)
       uvFade(:) = fadeInTable(inFadePos:endPos);
       synthFade(:) = fadeOutTable(inFadePos:endPos);
       outFadePos = endPos;
    else
        synthFade(:) = 0;
        uvFade(:) = 1;
        outFadePos = inFadePos;
    end
end
if (bvVUVbool && lvVUVbool)% if both BV and LV are voiced
           endPos = inFadePos-hopSize;

    if(endPos > 1)
       uvFade(:) = fadeInTable(inFadePos-1:-1:endPos);
       synthFade(:) = fadeOutTable(inFadePos-1:-1:endPos);
       outFadePos = endPos;
    else
        synthFade(:) = 1;
        uvFade(:) = 0;
        outFadePos = inFadePos;
    end
end    
end
    




