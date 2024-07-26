function [names,corr_out,p_out] = rmpWFCorr(filePath,saveFlag)
% LKW 10/19/21
% Inputs: 
%   filePath = string path to saved .mat file containing data
% Outputs:

corrCell = who('-file', filePath);
load(filePath);

nRamps = 21;
rampLim = 1;
rampPercs = linspace(0,rampLim,nRamps);
ct = 1;

for i = 1:size(corrCell,1)
    datVect = eval(corrCell{i,1});
    if isa(datVect,'double')
        tmpXs = linspace(0,1,length(datVect));
%         lenDatVect = length(datVect);     %In case datVect has missing extreme values
%         [tmpRho, tmpP] = corr([rampPercs(1:lenDatVect); datVect]');
        [tmpRho, tmpP] = corr([tmpXs; datVect]');
        corrCell(ct,2) = {tmpRho(1,2)};    %Get correlation value
        corrCell(ct,3) = {tmpP(1,2)};      %Get correlation p-val
        ct = ct + 1;
    elseif isa(datVect,'struct')
        tmpFields = fieldnames(datVect);
        for j = 1:size(tmpFields,1)
            tmpDat = datVect.(tmpFields{j});
            tmpXs = linspace(0,1,length(tmpDat));
            [tmpRho, tmpP] = corr([tmpXs; tmpDat]');
            corrCell(ct,1) = {tmpFields{j}};
            corrCell(ct,2) = {tmpRho(1,2)};    %Get correlation value
            corrCell(ct,3) = {tmpP(1,2)};      %Get correlation p-val
            ct = ct + 1;
        end
    end
end

names = corrCell(:,1);
corr_out = cell2mat(corrCell(:,2));
p_out = cell2mat(corrCell(:,3));

if saveFlag == 1
    sname = [filePath,'_corrs'];
    save(sname,'names','corr_out','p_out');
end

end