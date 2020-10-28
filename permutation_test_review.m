function[] = permutation_test_review(which_strain, which_nucleus)
%e.g.: permutation_test_review('RC',{'dLGN','vLGN','OPN','pret'})

rng(1); %for reproducibility
Nshuffle = 100000;
filepath = 'Data\classification\';
filename = []; count_nuclei = 0;  count = 0;
x_infra = []; x_gamma = []; x_both = []; x_mfr = []; x_rec = [];
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
        x_mfr = [x_mfr log10(mfr{m})];
        x_infra = [x_infra  is_infra{m}==1];
        x_gamma = [x_gamma  is_gamma{m}==1];
        x_both = [x_both  (is_infra{m} & is_gamma{m})];
        x_rec = [x_rec count*ones(1,numel(mfr{m}))];
    end
end
%divide into firing rate tiers
hfr = [-2:1:2]; 
[~,ind_mfr] = histc(x_mfr,hfr);
%statistic in the original dataset
stat = sum(x_both);
%create shuffle sets
stat_shuffle = zeros(1, Nshuffle);
for n = 1:Nshuffle
    stat_shuffle(n) = get_stat_shuffle(ind_mfr, x_rec, x_infra, x_gamma);
end
%compare stat and stat_shuffle get p-value
pval = 1-(sum(stat>stat_shuffle)/Nshuffle);
%make figure
fig = figure;
set(fig,'Position',[300 300 400 300]);
hxmin = min(min(stat_shuffle),stat);
hxmax = max(max(stat_shuffle),stat);
Nh = 12; 
hx = [hxmin:(hxmax-hxmin)/(Nh-1):hxmax];
hy = hist(stat_shuffle,hx);
subp = subplot(1,1,1); hold on;
bar(hx,hy,'FaceColor',[0.6 0.6 0.6]);
line([stat stat],[0 1.1*max(hy)],'Color','k','LineWidth',2);
title(['p-value = ' num2str(pval)]);
set(subp,'FontSize',12);
xlabel('#infra+beta/gamma','FontSize',14);
ylabel('#permutation','FontSize',14);
l = legend('shuffle','original');
set(l,'FontSize',12);

function[stat_shuffle] = get_stat_shuffle(mfr, rec, infra, hfo)
which_rec = unique(rec);
Nrec = numel(which_rec);
mfr_val = unique(mfr); Nmfr = numel(mfr_val);
stat_shuffle = 0; 
for n = 1:Nrec
    for m = 1:Nmfr
        ind = find((rec == which_rec(n))&(mfr==mfr_val(m)));
        Nind = numel(ind);
        infra(ind) = infra(ind(randperm(Nind))); 
        hfo(ind) = hfo(ind(randperm(Nind))); 
    end
end
stat_shuffle = sum(infra&hfo);


