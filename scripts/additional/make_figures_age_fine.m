
function make_figures_age_fine(InputData, NameVariable, TitlePlot, NameRepeated, Measure) 

close all

% InputData should be in format URSI, Visit, VariableOfInterest, Age

% script to plot age against FC measures (both raw values as well as
% residualized with respect to other terms in model)

WideFormat = unstack(InputData, NameVariable, NameRepeated); % this returns a file in wide with URSI, age and 3 columns for F9 per visit
WideFormat.MeanVariableOfInterest = nanmean(table2array(WideFormat(:,[3,4,5])),2); % compute mean across visits

MeanVariableOfInterest = table2array(WideFormat(:,6)); % extract meanF9 variable for plotting
Age_wide = table2array(WideFormat(:,2)); % extract age

figure; scatter(Age_wide, MeanVariableOfInterest)
lsline

% plot some
title(TitlePlot);
xlabel('Age (years/months)');
ylabel(Measure);
print(['Age_fine_effects ' TitlePlot], '-djpeg');

end