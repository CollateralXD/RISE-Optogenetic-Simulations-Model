function [handle] = plot_CA3LearnAct(phaseFlag,learnActivity,learnIn)

N = size(learnActivity,1);
handle = figure; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]); 

cPyr = cool(N);
cAIn = gray(N+2);
for i = 1:N
    plot(learnActivity(i,1:end-1), 'Color', cPyr(i,:)); hold on;
    plot(learnIn(i,:)*5, 'Color', cAIn(i,:));
end
xlabel('Time (ms)'); ylabel('Activation');
set(gca,'fontname','times','FontSize',24);
if phaseFlag == 1
    title('Learn Phase')
else
    title('Test Phase')

end