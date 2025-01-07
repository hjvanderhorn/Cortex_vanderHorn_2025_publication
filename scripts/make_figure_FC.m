function make_figure_FC(Data, group, visit, title_name, yaxis_name, PositionLegend, ImageFormat)

N_V1_HC = sum(group == 1 & visit == 1);
N_V2_HC = sum(group == 1 & visit == 2);
N_V3_HC = sum(group == 1 & visit == 3);
N_V1_TBI = sum(group == 2 & visit == 1);
N_V2_TBI = sum(group == 2 & visit == 2);
N_V3_TBI = sum(group == 2 & visit == 3);

means_Data_HC = [mean(Data(group == 1 & visit == 1)) ...
    mean(Data(group == 1 & visit == 2)) ...
    mean(Data(group == 1 & visit == 3))];

std_Data_HC = [std(Data(group == 1 & visit == 1)) ...
    std(Data(group == 1 & visit == 2)) ...
    std(Data(group == 1 & visit == 3))];

errors_Data_HC = std_Data_HC./ sqrt([N_V1_HC N_V2_HC N_V3_HC]); % get sem

means_Data_TBI= [mean(Data(group == 2 & visit == 1)) ...
    mean(Data(group == 2 & visit == 2)) ...
    mean(Data(group == 2 & visit == 3))];

std_Data_TBI = [std(Data(group == 2 & visit == 1)) ...
    std(Data(group == 2 & visit == 2)) ...
    std(Data(group == 2 & visit == 3))];

errors_Data_TBI = std_Data_TBI./ sqrt([N_V1_TBI N_V2_TBI N_V3_TBI]); % get sem


% errorbar(means_Data_HC, errors_Data_HC); hold on;
% errorbar(means_Data_TBI, errors_Data_TBI);

s=shadedErrorBar([1:3], means_Data_HC, errors_Data_HC, 'lineprops','-b');
set(s.edge,'LineWidth',1);
s.mainLine.LineWidth = 1.25;

hold on; 
s2=shadedErrorBar([1:3], means_Data_TBI, errors_Data_TBI, 'lineprops','-r');
set(s2.edge,'LineWidth',1);
s2.mainLine.LineWidth = 1.25;

%xlim([0 4]);
xticks([1:1:3]);
xticklabels({'~1 week', '~4 months', '~1 year'});

ax=gca;

% to make sure that there is a little bit of room before scale begins:
axis(ax, 'tight');
xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05);
ylim(ax, ylim(ax) + [-1,1]*range(ylim(ax)).* 0.05);

ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';

legend({'HC', 'pmTBI'}, 'Location', PositionLegend, 'FontSize',12);
title(title_name, 'FontSize', 14);
xlabel([]);
ylabel(yaxis_name, 'FontSize', 14, 'FontWeight', 'bold');

hold on;
% some additional plots to visualize data as points
plot([1:3], means_Data_HC, '.', 'Color', 'b', 'MarkerSize',20);
hold on;
plot([1:3], means_Data_TBI, '.', 'Color', 'r', 'MarkerSize',20);

ax.Legend.String(3:4) = []; % remove additional names from legend that are associated with additional plots

if ~strcmp(ImageFormat, 'none')
cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\figures');
print(title_name, ImageFormat,  '-r300');
end

end
