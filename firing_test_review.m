function[] = firing_test_review(which_strain, which_nucleus)
%e.g.: firing_test_review('RC',{'dLGN','vLGN','OPN','pret'})

filepath = 'Data\classification\';
filename = []; count_nuclei = 0;  count = 0;
x_infra = []; x_gamma = []; 
x_both = []; x_none = [];
x_mfr = [];
Nnuclei = numel(which_nucleus);
%
for n = 1:Nnuclei
    count_nuclei = count_nuclei+1;
    %load data
    filename_in = [which_nucleus{count_nuclei} '_' which_strain '_infra_calc_res']; 
    load([filepath filename_in],'is_infra','mfr');
    filename_gm = [which_nucleus{count_nuclei} '_' which_strain '_gamma_calc_res']; 
    load([filepath filename_gm],'is_gamma');
    %fill up xs
    for m = 1:numel(is_infra)
        count = count+1;
        x_none = [x_none ((~is_infra{m})&(~is_gamma{m}))];
        x_infra = [x_infra ((is_infra{m})&(~is_gamma{m}))];
        x_gamma = [x_gamma ((~is_infra{m})&is_gamma{m})];
        x_both = [x_both (is_infra{m}&is_gamma{m})];
        x_mfr = [x_mfr mfr{m}];
    end
end
fr_none = (x_mfr(find(x_none))); Nnone = sum(x_none);
fr_infra = (x_mfr(find(x_infra))); Ninfra = sum(x_infra);
fr_gamma = (x_mfr(find(x_gamma))); Ngamma = sum(x_gamma);
fr_both = (x_mfr(find(x_both))); Nboth = sum(x_both);
hfr = (-1:0.5:2);

%
% figure; hold on;
% plot(hfr,hist(log10(fr_none),hfr)/Nnone,'k','LineWidth',2);
% plot(hfr,hist(log10(fr_infra),hfr)/Ninfra,'b','LineWidth',2);
% plot(hfr,hist(log10(fr_gamma),hfr)/Ngamma,'r','LineWidth',2);
% plot(hfr,hist(log10(fr_both),hfr)/Nboth,'g','LineWidth',2);

fig = figure; 
set(fig,'Position',[100 100 400 350]);
hold on;
p_none = hist(log10(fr_none),hfr)/Nnone;
p_infra = hist(log10(fr_infra),hfr)/Ninfra;
p_gamma = hist(log10(fr_gamma),hfr)/Ngamma;
p_both =  hist(log10(fr_both),hfr)/Nboth;
%
plot(hfr,p_none,'k.-','LineWidth',2,'MarkerSize',24);
plot(hfr,p_infra,'b.-','LineWidth',2,'MarkerSize',24);
plot(hfr,p_gamma,'r.-','LineWidth',2,'MarkerSize',24);
plot(hfr,p_both,'g.-','LineWidth',2,'MarkerSize',24);
legend({'none','infra-only','gamma-only','infra&gamma'});
xlabel('log10(FR(Hz))'); ylabel('Probability');
%
fig = figure; 
set(fig,'Position',[100 100 400 350]);
h = subplot(1,1,1); hold on;
M = [mean(fr_none) mean(fr_infra) mean(fr_gamma) mean(fr_both)];
%S = [std(fr_none)/sqrt(Nnone) std(fr_infra)/sqrt(Ninfra) std(fr_gamma)/sqrt(Ngamma) std(fr_both)/sqrt(Nboth)];
S = [std(fr_none) std(fr_infra) std(fr_gamma) std(fr_both)];
bar(M,'BarWidth',0.5,'FaceColor',0.666*ones(1,3));
errorbar(M,S,'.k','LineWidth',2);
set(h,'XTick',1:4);
set(h,'XTickLabel',{'-','I','G','I&G'});
p = kruskalwallis([fr_none fr_infra fr_gamma fr_both],[1*ones(1,Nnone) 2*ones(1,Ninfra) 3*ones(1,Ngamma) 4*ones(1,Nboth)],'off');
ylabel('FR(Hz)');
title(p);
%compare individual groups
p_none_infra = ranksum(fr_none,fr_infra);
p_infra_gamma = ranksum(fr_infra,fr_gamma);
p_gamma_both = ranksum(fr_gamma,fr_both);
disp(sprintf('p none-infra: %s ; p infra-gamma: %s ; p gamma-both %s',num2str(p_none_infra),num2str(p_infra_gamma),num2str(p_gamma_both))) 