function[] = infra_phase_gamma_pow_review(which_genotype)
%e.g. infra_phase_gamma_pow_review('RC')

phase_range = [-pi:pi/2:pi];
count_same = 0; 
count_distinct = 0;
%
count = 1;
filename = [];
if strcmp(which_genotype,'RC')
    filename{count} = 'dLGN_RC'; count = count+1;
    filename{count} = 'vLGN_RC'; count = count+1;
    filename{count} = 'OPN_RC'; count = count+1;
    filename{count} = 'pret_RC'; count = count+1;
elseif strcmp(which_genotype,'MELKO')
    filename{count} = 'dLGN_MELKO'; count = count+1;
    filename{count} = 'vLGN_MELKO'; count = count+1;
    filename{count} = 'OPN_MELKO'; count = count+1;
    filename{count} = 'pret_MELKO'; count = count+1;
else
    return;
end
%
Nfile = numel(filename); 
filepath1 = 'Data\';
filepath2 = 'Data\classification\';
%
isinfra = []; xb = []; isgamma = [];
freq_infra = []; freq_gamma = [];
for n = 1:Nfile
    load([filepath1 filename{n} '_data'],'x');
    xb = [xb x];
    load([filepath2 filename{n} '_infra_calc_res'],'is_infra','freq');
    isinfra = [isinfra is_infra];
    freq_infra = [freq_infra freq];
    load([filepath2 filename{n} '_gamma_calc_res'],'is_gamma','freq');
    isgamma = [isgamma is_gamma]; 
    freq_gamma = [freq_gamma freq];
end
N = numel(isinfra);
%
for n = 1:numel(xb)
    %define units with infra and/or gamma
    temp_infra = isinfra{n}; Finfra = freq_infra{n};
    temp_gamma = isgamma{n}; Fgamma = freq_gamma{n};
    ind_infra_only = find(temp_infra & ~temp_gamma);
    ind_gamma_only = find(~temp_infra & temp_gamma);
    ind_infra_gamma = find(temp_infra & temp_gamma);
    Ninfra_only = numel(ind_infra_only);
    Ngamma_only = numel(ind_gamma_only);
    Ninfra_gamma = numel(ind_infra_gamma);
    %look at uncoupled PAC
    if((Ninfra_only>0)&(Ngamma_only>0))
        %pool infra
        xphase = pool_bin_filter_transform(xb{n},ind_infra_only,'infra',Finfra(ind_infra_only));
        %pool gamma
        xabs = pool_bin_filter_transform(xb{n},ind_gamma_only,'gamma',Fgamma(ind_gamma_only));
        %resample gamma
        xabs = resample(xabs,numel(xphase),numel(xabs));
        %calculate mean xabs as function of phase
        count_distinct = count_distinct+1;
        xabs_phase_distinct(count_distinct,:) = calculate_amplitude_phase(xabs,xphase,phase_range);
    end
    if(Ninfra_gamma>0)
        %pool infra
        xphase = pool_bin_filter_transform(xb{n},ind_infra_gamma,'infra',Finfra(ind_infra_gamma));
        %pool gamma
        xabs = pool_bin_filter_transform(xb{n},ind_infra_gamma,'gamma',Fgamma(ind_infra_gamma));
        %resample gamma
        xabs = resample(xabs,numel(xphase),numel(xabs));
        %calculate mean xabs as function of phase
        count_same = count_same+1;
        xabs_phase_same(count_same,:) = calculate_amplitude_phase(xabs,xphase,phase_range);
    end
    disp(sprintf('rec %s of %s',num2str(n),num2str(numel(xb))));
end
%
stat = mean(xabs_phase_same(:,[2 3]),2) - mean(xabs_phase_same(:,[1 4]),2);
psame = signtest(stat);
stat = mean(xabs_phase_distinct(:,[2 3]),2) - mean(xabs_phase_distinct(:,[1 4]),2);
pdistinct = signtest(stat);
%
save([filepath1 which_genotype '_PAC_results'],'xabs_phase_same','xabs_phase_distinct','psame','pdistinct');
%
figure; hold on;
bar(mean(xabs_phase_same),'BarWidth',0.5,'FaceColor',0.666*ones(1,3)); 
errorbar(mean(xabs_phase_same),std(xabs_phase_same),'k.','LineWidth',2);
title(['same: pval=' num2str(psame)]);
%
figure; hold on;
bar(mean(xabs_phase_distinct),'BarWidth',0.5,'FaceColor',0.666*ones(1,3)); 
errorbar(mean(xabs_phase_distinct),std(xabs_phase_distinct),'k.','LineWidth',2);
title(['distinct: pval=' num2str(pdistinct)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[xout] = pool_bin_filter_transform(x,ind,which_band,freq)
freq = mean(freq);
if strcmp(which_band,'infra')
    dT = 1; T = [0:dT:900];
    c = kaiserord([2.5*freq 5*freq], [1 0], [0.025 0.2], 1/dT, 'cell');
    b = fir1(c{:});
elseif strcmp(which_band,'gamma')
    dT = 0.0025; T = [0:dT:900];
    c = kaiserord([0.6*freq  0.8*freq 1.2*freq 1.4*freq], [0 1 0], [0.1 0.01 0.1], 1/dT, 'cell');
    b = fir1(c{:});
else return
end
Nt = numel(T);
Nind = numel(ind); 
xbin = zeros(1,Nt);
for n = 1:Nind
    temp = x(:,ind(n));
    temp = temp(temp>0);
    xbin = xbin+hist(temp,T);
end
%detrend
xbin = detrend(xbin);
%pass-band filter
xbin_filt = filtfilt(b,1,xbin);
if strcmp(which_band,'infra')
    %get phase
    xout = wrapToPi(phase(hilbert(xbin_filt)));
elseif strcmp(which_band,'gamma') 
    %get amplitude
    xout = abs(hilbert(xbin_filt));
    %low pass amplitude
    c = kaiserord([0.05 1], [1 0], [0.01 0.1], 200, 'cell');
    b = fir1(c{:});
    xout = filtfilt(b,1,xout);
    xout = xout(1:20:end);
else return
end


function[xabs_phase] = calculate_amplitude_phase(xabs,xphase,phase_range)
%reduce
xabs = xabs(25:end-24);
xphase = xphase(25:end-24);
%bin phase values 
[~,phase_val] = histc(xphase,phase_range); 
Nphase = numel(phase_range)-1;
%calculate norm amplitude as function of phase
xabs_phase = zeros(1,Nphase);
for n = 1:Nphase
    xabs_phase(n) = mean(xabs(phase_val==n));
end
xabs_phase = xabs_phase/sum(xabs_phase);
