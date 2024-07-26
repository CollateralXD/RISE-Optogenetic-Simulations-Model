function [dOut,meansOut,sdsOut] = learn_CohenStats(actCell,hactCell,pStruct)
%Run Stats on pyramidals only for the CA3Net Ripple Extend code
%Input: actCell = cell containing the activity of pyramidals for ramp,
%square and control
%       hactCell = cell with activity of interneurons
% 
%Output: cohen's d of the comparison to control and raw mean and stdev of
%opto ramp
%LKW 4/24/21
%%%

aRamp = actCell{1}; aControl = actCell{2};
% hRamp = hactCell{1}; hControl = hactCell{2};

for  i = 1:pStruct.N    %Find the highest peak, or return no peak if flat line
    [pksTmp,locsTmp] = findpeaks(aRamp(i,:));
    if ~isempty(locsTmp); pksR(i) = max(pksTmp); locsR(i) = locsTmp(pksTmp == pksR(i)); end
    [pksTmp,locsTmp] = findpeaks(aControl(i,:));
    if ~isempty(locsTmp); pksC(i) = max(pksTmp); locsC(i) = locsTmp(pksTmp == pksC(i)); end
end

locsCell = {locsR,locsC};
pksCell = {pksR,pksC};

%Threshold method
actThresh = pStruct.actThresh;
tttCell = {};

for i = 1:2
    pksTmp = pksCell{i};        %All the peaks for that sim (e.g. ramp)
    actTmp = actCell{i};        %Pyramidal activation for that sim e.g. NxT
    ttt = [];                   %Vector of threshold timestamps
    for j = 1:numel(pksTmp)     %For each peak in the sim
        indtmp = find(actTmp(j,:)>actThresh(j),1);  %Get index of threshold crossing point
        if ~isempty(indtmp)
            ttt = [ttt,indtmp]; %If index exists add that timestamp to vector
        end
    end
    tttCell(i) = {ttt};         %Cell of vectors for timestamps of threshold
end

%Prep for stats
rampDeltas = diff(tttCell{1});
ctlDeltas = diff(tttCell{2});

meansOut = mean(rampDeltas);
sdsOut = std(rampDeltas);

try
    dOut = computeCohen_d(rampDeltas,ctlDeltas,'independent');
catch 
    dOut = NaN;
end

end