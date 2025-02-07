function [] = post_silhouette_plot(theta_samples,l_post_samples, par_names, par_num)
%The function "post_silhouette_plot" displays the silhouette of the 
%posterior surface from the chosen parameter's perspective

par = theta_samples(par_num,:);

post_samples = exp(l_post_samples/10);
%post_samples = (l_post_samples)/2;

[par_sort, ind] = sort(par);
post_samples_sort = post_samples(ind);

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot(par_sort, post_samples_sort)

xlabel(par_names(1,par_num), 'FontSize', 12), ylabel('Posterior Density', 'FontSize', 12)
fprintf("\nmaxval = %d",max(post_samples_sort)*1.1);
fprintf("\nminval = %d",min(post_samples_sort)*1.1);
% ylim([min(0,max(post_samples_sort)*1.1) max(0,max(post_samples_sort)*1.1)])

hold off

end