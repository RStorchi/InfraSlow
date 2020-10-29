function[] = FF_figure_review()
filepath1 = 'Data\';
load([filepath1 'FF_calc'],'T','FF','infra');
%
FFinfra = FF(infra==1,:); 
FFinfra_min = quantile(FFinfra,0.01);
FFinfra_max = quantile(FFinfra,0.99);
FFnot = FF(infra==0,:);
FFnot_min = quantile(FFnot,0.01);
FFnot_max = quantile(FFnot,0.99);
%
figure; hold on;
line([0.001 100],[1 1],'Color',0.666*ones(1,3),'LineWidth',2);
errorbar(T,mean(FFnot), FFnot_min, FFnot_max,'k','LineWidth',2);
errorbar(T,mean(FFinfra), FFinfra_min, FFinfra_max,'b','LineWidth',2);

