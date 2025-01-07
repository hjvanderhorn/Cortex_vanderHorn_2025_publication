
function make_figures_age(InputData, NameVariable, TitlePlot, NameRepeated, Measure, ImageFormat) 


% InputData should be in format URSI, Visit, VariableOfInterest, Age

% script to plot age against FC measures (both raw values as well as
% residualized with respect to other terms in model)

WideFormat = unstack(InputData, NameVariable, NameRepeated); % this returns a file in wide with URSI, age and 3 columns for F9 per visit
WideFormat.MeanVariableOfInterest = nanmean(table2array(WideFormat(:,[3,4,5])),2); % compute mean across visits

MeanVariableOfInterest = table2array(WideFormat(:,6)); % extract meanF9 variable for plotting
Age_wide = table2array(WideFormat(:,2)); % extract age

AgeGroups = unique(Age_wide); % make age groups

% now compute the mean and standard error per age group
for i=1:size(AgeGroups,1) 
    M(i) = mean(MeanVariableOfInterest(Age_wide==AgeGroups(i))); 
    SEM(i) = std(MeanVariableOfInterest(Age_wide==AgeGroups(i)))./sqrt(sum(Age_wide==AgeGroups(i))); % sem
end

% third order polynomial fit:
p = polyfit(AgeGroups', M,3);
y = polyval(p, AgeGroups);

% plot some
figure; %errorbar(AgeGroups,M, SE);
s=shadedErrorBar(AgeGroups,M, SEM);

set(s.edge,'LineWidth',1);
s.mainLine.LineWidth = 1.25;
s.mainLine.LineStyle = 'none';
hold on;
plot(AgeGroups, y,'k', 'LineWidth', 1.25);
hold on; plot(AgeGroups, M, ':.', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 0.01 );

title(TitlePlot, 'FontSize', 14);
xlabel('Age (years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(Measure, 'FontSize', 14, 'FontWeight', 'bold');
xticks(min(AgeGroups):1:max(AgeGroups));

ax=gca;
axis(ax, 'tight');
xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05);
ylim(ax, ylim(ax) + [-1,1]*range(ylim(ax)).* 0.05);
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

if ~strcmp(ImageFormat, 'none')
cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\figures');
print(['Age effects ' TitlePlot],ImageFormat,  '-r300' );
end

% 
% figure; %errorbar(AgeGroups,smoothdata(M), SE);
% s2=shadedErrorBar(AgeGroups,smoothdata(M), SE);
% set(s2.edge,'LineWidth',1);
% s2.mainLine.LineWidth = 1.25;
% title([TitlePlot ', smoothed curve'], 'FontSize', 16);
% xlabel('Age (years)', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel(Measure, 'FontSize', 14, 'FontWeight', 'bold');
% xticks(min(AgeGroups):1:max(AgeGroups));
% 
% ax=gca;
% axis(ax, 'tight');
% xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05);
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% 
% print(['Age effects ' TitlePlot ', smoothed curve'], '-djpeg');

end