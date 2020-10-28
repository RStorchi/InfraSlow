function[] = gamma_classification_preprocess(filename)
%load data
filepath = 'Data\';
load([filepath filename '_data']);
%bin data
dT = 0.005; fs = 1/dT; pow = 9;
wind = 2^pow; noverlap = 2^(pow-1); nfft = 2^pow;
Tbin = 0:dT:900; Nt = numel(Tbin);
xbin = cell(1,Nsheet);
xpsd =  cell(1,Nsheet);
xmfr = cell(1,Nsheet);
xok = cell(1,Nsheet);
N = 0;
for n = 1:Nsheet
    for m = 1:size(x{n},2)
        xbin{n}(:,m) = hist(x{n}(:,m),Tbin); 
        [xpsd{n}(:,m),F] = pwelch(detrend(xbin{n}(:,m)),wind,noverlap,nfft,fs);
        xpsd{n}(:,m) = xpsd{n}(:,m)/sum(xpsd{n}(:,m)); %don't normalize here
        xmfr{n}(m) = mean(xbin{n}(:,m))/dT;
        xok{n}(m) = ~isnan(xpsd{n}(1,m));
        N = N+1;
        disp(sprintf('cell %s',num2str(N)));
    end
end
%save
save([filepath filename '_gamma_calc'],'xpsd','xmfr','xok','F');

