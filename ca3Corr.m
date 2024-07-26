function [corrCell] = ca3Corr(filePath)
% LKW 4/2/21
% Inputs: 
%   xVect = 1D vector of independent variable, probably rampPercs
%   DatVect = 1D vector of data to get correlations
% Outputs:

corrCell = who('-file', filePath);
load(filePath);

nRamps = 21;
rampLim = 1;
rampPercs = linspace(0,rampLim,nRamps);
ct = 1;

for i = 1:size(corrCell,1)
    datVect = eval(corrCell{i,1});
    lenDatVect = length(datVect);     %In case datVect has missing extreme values
%     if isstruct(datVect)
%         strucCell = fieldnames(datVect);
%         for j = 1:size(strucCell,1)
%             [tmpRho, tmpP] = corr([rampPercs(1:lenDatVect); datVect]');
%             corrCell(ct,2) = {tmpRho(1,2)};    %Get correlation value
%             corrCell(ct,3) = {tmpP(1,2)};      %Get correlation p-val
%             ct = ct + 1;
%         end
%     else
        [tmpRho, tmpP] = corr([rampPercs(1:lenDatVect); datVect]');
        corrCell(ct,2) = {tmpRho(1,2)};    %Get correlation value
        corrCell(ct,3) = {tmpP(1,2)};      %Get correlation p-val
        ct = ct + 1;
%     end
end

end