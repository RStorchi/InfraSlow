function[] = infra_crosscorr_review(which_nucleus, which_strain)
%e.g. infra_crosscorr_review('pret','MELKO')

%load data
filepath1 = 'Data\';
filepath2 = 'Data\classification\';
%load  spike counts
filename_in = [which_nucleus '_' which_strain '_infra_calc']; 
load([filepath1 filename_in],'xbin');
%load classifications
filename_in = [which_nucleus '_' which_strain '_infra_calc_res']; 
load([filepath2 filename_in],'is_infra');
filename_in = [which_nucleus '_' which_strain '_gamma_calc_res']; 
load([filepath2 filename_in],'is_gamma');

%get infra_gamma indexes
Nrec = numel(xbin);
is_infra_gamma = cell(1,Nrec); 
for n = 1:Nrec
    is_infra_gamma{n} = is_infra{n} & is_gamma{n};
    is_infra_only{n} = is_infra{n} & (~is_gamma{n});
end

%calculate x-corr
max_lag = 300; 
xc_infra_only = cell(1,Nrec); count_infra_only = 0; 
xc_infra_gamma = cell(1,Nrec); count_infra_gamma = 0; 
xc_between = cell(1,Nrec); count_between = 0; 
for n = 1:Nrec
    xbin{n} = detrend(xbin{n}); 
    %xc among infra-only
    ind_infra_only = find(is_infra_only{n});
    if numel(ind_infra_only)>1
        [xc_infra_only{n},tc] = get_infra_crosscorr_review_v1_sub1(xbin{n}(:,ind_infra_only),max_lag);
    end
    %xc among infra&gamma
    ind_infra_gamma = find(is_infra_gamma{n});
    if numel(ind_infra_gamma)>1
        xc_infra_gamma{n} = get_infra_crosscorr_review_v1_sub1(xbin{n}(:,ind_infra_gamma),max_lag);
    end
    %xc between infra-only and infra&gamma
    if ((numel(ind_infra_gamma)>1) & (numel(ind_infra_gamma)>1))
        xc_between{n} = get_infra_crosscorr_review_v1_sub2(xbin{n}(:,ind_infra_only),xbin{n}(:,ind_infra_gamma),max_lag);
    end
end
%merge recordings
xc_infra_only = vertcat(xc_infra_only{:}); Ninfra_only = size(xc_infra_only,1);
xc_infra_gamma = vertcat(xc_infra_gamma{:}); Ninfra_gamma = size(xc_infra_gamma,1);
xc_between = vertcat(xc_between{:}); Nbetween = size(xc_between,1);
%stats
indt = find((tc>=-10)&(tc<=10));
if Ninfra_only
    [mxc_infra_only,pval_infra_only] = get_infra_crosscorr_review_v1_sub3(xc_infra_only,indt);
    disp(sprintf('infra_only: r=%s pval=%s n=%s',num2str(mxc_infra_only),num2str(pval_infra_only),num2str(Ninfra_only))); 
end
if Ninfra_gamma
    [mxc_infra_gamma,pval_infra_gamma] = get_infra_crosscorr_review_v1_sub3(xc_infra_gamma,indt);
    disp(sprintf('infra_gamma: r=%s pval=%s n=%s',num2str(mxc_infra_gamma),num2str(pval_infra_gamma),num2str(Ninfra_gamma))); 
end
if Nbetween
    [mxc_between,pval_between] = get_infra_crosscorr_review_v1_sub3(xc_between,indt);
    disp(sprintf('between: r=%s pval=%s n=%s',num2str(mxc_between),num2str(pval_between),num2str(Nbetween)));
end
%
figure; hold on;
plot(tc,mean(xc_infra_only),'LineWidth',2);
plot(tc,mean(xc_infra_gamma),'LineWidth',2);
plot(tc,mean(xc_between),'LineWidth',2);
legend({'infra-only','infra&gamma','between'});

function[xc,tc] = get_infra_crosscorr_review_v1_sub1(x,max_lag)
xc = []; count = 0; 
N = size(x,2); 
for n = 1:N
    for m = n+1:N
        count = count+1;
        [xc(count,:),tc] = xcorr(x(:,n),x(:,m),'coeff',max_lag);
    end
end
    
function[xc,tc] = get_infra_crosscorr_review_v1_sub2(x,y,max_lag)
xc = []; count = 0; 
Nx = size(x,2); Ny = size(y,2); 
for n = 1:Nx
    for m = 1:Ny
        count = count+1;
        [xc(count,:),tc] = xcorr(x(:,n),y(:,m),'coeff',max_lag);
    end
end

function[mmxc,pval] = get_infra_crosscorr_review_v1_sub3(xc,indt)
N = size(xc,1);
mxc = median(xc(:,indt)');
mmxc = mean(mxc);
pval = signtest(mxc);

