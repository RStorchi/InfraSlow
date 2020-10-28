function[]=gamma_classification_review(is_graph)

count = 1;
filename = [];
filename{count} = 'dLGN_RC_gamma_calc'; count = count+1;
filename{count} = 'dLGN_MELKO_gamma_calc'; count = count+1;
filename{count} = 'dLGN_RDCL_gamma_calc'; count = count+1;

filename{count} = 'vLGN_RC_gamma_calc'; count = count+1;
filename{count} = 'vLGN_MELKO_gamma_calc'; count = count+1;
filename{count} = 'vLGN_RDCL_gamma_calc'; count = count+1;

filename{count} = 'OPN_RC_gamma_calc'; count = count+1;
filename{count} = 'OPN_MELKO_gamma_calc'; count = count+1;
filename{count} = 'OPN_RDCL_gamma_calc'; count = count+1;

filename{count} = 'pret_RC_gamma_calc'; count = count+1;
filename{count} = 'pret_MELKO_gamma_calc'; count = count+1;
filename{count} = 'pret_RDCL_gamma_calc'; count = count+1;
%
Nfile = numel(filename); 
filepath = 'Data\';
%
psd_norm = []; mfr = []; is_ok = [];
Nsheet = zeros(1,Nfile); Nunit = cell(1,Nfile);
for n = 1:Nfile
    load([filepath filename{n}],'xpsd','xmfr','xok','F');
    psd_norm = [psd_norm horzcat(xpsd{:})];
    mfr = [mfr horzcat(xmfr{:})];
    is_ok = [is_ok horzcat(xok{:})];
    Nsheet(n) = numel(xmfr);
    for m = 1:Nsheet(n)
        Nunit{n} = [Nunit{n} size(xmfr{m},2)];
    end
end
N = numel(mfr);
%%%%%%%%%%%%%%%%%%%%%%%gamma criterion%%%%%%%%%%%%%%%%%%%%%%%
is_gamma = false(1,N);
Fmax = zeros(1,N);
indf = find((F>20)&(F<100));
TH = 0.;
count = 0;

for n = 1:N
    %psd_norm(:,n) = detrend(psd_norm(:,n));
    temp = (psd_norm(:,n)); 
    %temp_lf = filtfilt(ones(1,5),1,psd_norm(:,n));
    temp_lf = medfilt1(temp,20);
    temp0 = temp-temp_lf;
    temp =  temp0(indf);
    %THvar = quantile(abs(diff(temp)),0.95);
    %THvar = 6*std(abs(diff(temp)));
    TH = 5*std(temp);
    %temp = (psd_norm(indf,n)); 
    %[peaks, ind_peaks] = findpeaks(temp,indf,'Threshold',0,'SortStr','descend','MinPeakProminence',THvar);
    [peaks, ind_peaks] = findpeaks(temp,indf,'Threshold',0,'SortStr','descend');
    if numel(peaks)
        Fmax(n) = F(ind_peaks(1));
        if peaks(1)>TH
            is_gamma(n) = true;
        end
    end
    disp(sprintf('cell %s out of %s',num2str(n),num2str(N)));
    if is_graph & is_gamma(n)
        hold on;
        plot(F,psd_norm(:,n)','LineWidth',2);
        plot(F,temp_lf,'LineWidth',2);
        plot(F,temp0,'LineWidth',2);
        line([F(1) F(end)],[TH TH],'Color','k','LineWidth',2);
        %line([F(1) F(end)],[THvar THvar],'Color',0.666*ones(1,3),'LineWidth',2);
        if numel(peaks)
            plot(Fmax(n),temp(ind_peaks(1)-indf(1)+1),'r.','MArkerSize',15);
        end
        title(num2str(is_gamma(n)));
        ginput(); clf;
    end
end
Ngamma = sum(is_gamma);
%%%%%%%figure%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
%
temp1 = psd_norm(:,is_gamma)';
[~,ind_sort] = sort(Fmax(is_gamma));
temp1 = temp1(ind_sort,:);
%
temp2 = psd_norm(:,~is_gamma)';
[~,ind_sort] = sort(Fmax(~is_gamma));
temp2 = temp2(ind_sort,:);
%
temp = [temp1; temp2];
p = pcolor(F,1:N,temp); set(p,'LineStyle','none')
line([F(1) F(end)],[1 1]*Ngamma,'LineWidth',2,'Color','r');
xlim([F(1) F(end)]); ylim([1 N]);
caxis([0.00 0.0075]);
title(Ngamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write the results files (mat)
%close all;
is_gamma_all = is_gamma; 
mfr_all = mfr;
count0 = 1; count1 = 0;
for n = 1:Nfile
    is_gamma = cell(1,Nsheet(n));
    freq = cell(1,Nsheet(n));
    mfr = cell(1,Nsheet(n));
    for m = 1:Nsheet(n)
        count1 = count1+Nunit{n}(m);
        ind = [count0:count1];
        is_gamma{m} = is_gamma_all([count0:count1]);
        freq{m} = Fmax(count0:count1);
        mfr{m} = mfr_all(count0:count1);
        count0 = count1+1;
    end
    save([filepath 'classification\' filename{n} '_res'],'is_gamma','freq','mfr');
end

