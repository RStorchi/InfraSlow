function[]= infra_classification_review(is_graph)
THerr = 0.2250; THbr = 0.6;
%
count = 1;
filename = [];
filename{count} = 'dLGN_RC_infra_calc'; count = count+1;
filename{count} = 'dLGN_MELKO_infra_calc'; count = count+1;
filename{count} = 'dLGN_RDCL_infra_calc'; count = count+1;

filename{count} = 'vLGN_RC_infra_calc'; count = count+1;
filename{count} = 'vLGN_MELKO_infra_calc'; count = count+1;
filename{count} = 'vLGN_RDCL_infra_calc'; count = count+1;

filename{count} = 'OPN_RC_infra_calc'; count = count+1;
filename{count} = 'OPN_MELKO_infra_calc'; count = count+1;
filename{count} = 'OPN_RDCL_infra_calc'; count = count+1;

filename{count} = 'pret_RC_infra_calc'; count = count+1;
filename{count} = 'pret_MELKO_infra_calc'; count = count+1;
filename{count} = 'pret_RDCL_infra_calc'; count = count+1;
%
Nfile = numel(filename); 
filepath = 'Data\';
%
b = []; infra = [];
err0 = []; err1 = []; err2 = [];
ac = []; ac_fit = [];
xb = []; Nunit = cell(1,Nfile);
mfr = [];
for n = 1:Nfile
    load([filepath filename{n}],'ball2','is_infra','errall1','errall2','xac','xac_fit','tac','xbin');
    b = [b horzcat(ball2{:})];
    err1 = [err1 horzcat(errall1{:})];
    err2 = [err2 horzcat(errall2{:})];
    ac = [ac; horzcat(xac{:})'];
    ac_fit = [ac_fit; horzcat(xac_fit{:})'];
    xb = [xb; horzcat(xbin{:})'];
    mfr = [mfr mean(horzcat(xbin{:}))];
    Nsheet(n) = numel(ball2);
    for m = 1:Nsheet(n)
        Nunit{n} = [Nunit{n} size(ball2{m},2)];
    end
end
N = size(b,2);

%
fig = figure; hold on;
set(fig,'Position',[200 200 400 400]);
b1 = [100 300 500 100 300 500 100 300 500; [100 300 500]*1.2 [100 300 500]*0.6 [100 300 500]*0.3];
lg = {'class boundary','infra','on boundary','null'};
plot(0.1:0.1:900,log10(0.6*(0.1:0.1:900)),'k','LineWidth',2);
plot(b1(1,[1:3]),log10(b1(2,[1:3])),'.','MarkerSize',30,'Color','b'); 
plot(b1(1,[4:6]),log10(b1(2,[4:6])),'.','MarkerSize',30,'Color','k');
plot(b1(1,[7:9]),log10(b1(2,[7:9])),'.','MarkerSize',30,'Color',0.666*ones(1,3));
tac1 = tac(tac>=0);
xlim([0 600]); ylim([1 3.5]); 
xlabel('Period(s)','FontSize',12);
ylabel('Decay(log10(s))','FontSize',12);
legend(lg);
fig = figure;
set(fig,'Position',[200 200 400 400]);
for n = 1:9
    subplot(3,3,n);
    if n < 4
        plot(tac1,cos(2*pi*tac1/b1(1,n)).*exp(-tac1/b1(2,n)),'LineWidth',2,'Color','b');
    elseif n < 7
        plot(tac1,cos(2*pi*tac1/b1(1,n)).*exp(-tac1/b1(2,n)),'LineWidth',2,'Color','k');
    else
        plot(tac1,cos(2*pi*tac1/b1(1,n)).*exp(-tac1/b1(2,n)),'LineWidth',2,'Color',0.666*ones(1,3));
    end
    xlim([0 900]);
end
%
b0 = b;
b1 = zeros(3,N);
x1 = b(1,:).*b(3,:); x2 = b(4,:).*b(6,:);
bratio = (x1-x2)./(x1+x2);
for n = 1:N
    if bratio(n) > 0
        b1(:,n) = b(1:3,n);
    else
        b1(:,n) = b(4:6,n);
        b0(:,n) = [b(4:6,n);b(1:3,n)];
    end
end
b = b1;
c = (err1-err2)./(err1);
bratio = abs(bratio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%criterion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_infra = find((b(2,:)>50)&(b(2,:)<450)&((b(3,:)./b(2,:))>THbr)&(c>THerr));
ind_none = setdiff(1:N,ind_infra); Nnone = numel(ind_none);
Ninfra = numel(ind_infra);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
if is_graph;
    fig = figure; set(fig,'Position',[100 100 400 400]);
    for n = 1:Ninfra
        ind = ind_infra(n);
        t = tac(tac>0);
        subplot(2,1,1); plot(xb(ind,:),'k','LineWidth',2); xlim([0 900]);
        subplot(2,1,2); hold on;
        plot(t,ac(ind,tac>0),'k',t,ac_fit(ind,:),'b','LineWidth',2); 
        plot(t,b0(1,ind)*cos(2*pi*t/b0(2,ind)).*exp(-t/b0(3,ind)),'r','LineWidth',2); 
        plot(t,b0(4,ind)*cos(2*pi*t/b0(5,ind)).*exp(-t/b0(6,ind)),'g','LineWidth',2);
        xlim([0 900]);
        title(sprintf('%s out of %s',num2str(n),num2str(Ninfra)));
        ginput(); clf;
    end
end
%
figure; hold on;
plot(b(2,ind_none),log10(b(3,ind_none)),'.','MarkerSize',12,'Color',0.666*ones(1,3));
plot(b(2,ind_infra),log10(b(3,ind_infra)),'.','MarkerSize',12,'Color','b');
plot(1:900,log10((1:900)*0.5),'k','LineWidth',2);
ylim([0 4]);
xlim([0 900]);
title(sprintf('infra %s of %s (%s perc)',num2str(Ninfra),num2str(N),num2str(100*Ninfra/N)));
xlabel('Period(s)','FontSize',12);
ylabel('Decay(log10(s))','FontSize',12);
%
figure; 
h = subplot(1,1,1); hold on;
err = (err1-err2)./err1;
M = [mean(err(ind_infra)) mean(err(ind_none))];
%S = [std(err(ind_infra))/sqrt(Ninfra) std(err(ind_none))/sqrt(Nnone)];
S = [std(err(ind_infra)) std(err(ind_none))];
bar([1 2],M,'BarWidth',0.5,'FaceColor',0.666*ones(1,3));
errorbar([1 2],M,S,'k.','LineWidth',2);
set(h,'XTick',[1 2],'XTickLabel',{'Infra','Null'},'FontSize',14);
ylabel('(ErrorNull-ErrorOsc)/(ErrNull)','FontSize',14);
xlim([0 3]);
title(ranksum(err(ind_infra),err(ind_none)));
%
figure; 
h = subplot(1,1,1); hold on;
br = b0(3,:)./b0(2,:); 
xtemp = linspace(-2,3,50);
temp = hist(log10(br),xtemp);
bar(xtemp,temp,'BarWidth',0.8,'FaceColor',0.666*ones(1,3));
line(log10(THbr)*[1 1],1.1*[0 max(temp)],'Color','k','LineWidth',2);
xlabel('bratio','FontSize',14); ylabel('#unit','FontSize',14);
xlim([-2 3]);
%
figure; 
h = subplot(1,1,1); hold on;
xtemp = linspace(-0.1,1,50);
temp = hist(c,xtemp);
bar(xtemp,temp,'BarWidth',0.8,'FaceColor',0.666*ones(1,3));
line(THerr*[1 1],1.1*[0 max(temp)],'Color','k','LineWidth',2);
xlabel('error','FontSize',14); ylabel('#unit','FontSize',14);
xlim([-0.1 1]);
%
figure;
h = subplot(1,1,1); hold on;
xtemp = linspace(-3,3,100);
temp = hist(log10(b0(2,:)./b0(5,:)),xtemp);
bar(xtemp,temp,'BarWidth',0.8,'FaceColor',0.666*ones(1,3));
line(log10(0.5)*[1 1],1.1*[0 max(temp)],'Color','k','LineWidth',2);
line(log10(2)*[1 1],1.1*[0 max(temp)],'Color','k','LineWidth',2);
xlabel('log10 freq ratio','FontSize',14); ylabel('#unit','FontSize',14);
xlim([-2 3]);
%
[~,ind_sort] = sort(b(2,ind_infra));
ac(ind_infra,:) = ac(ind_infra(ind_sort),:);
[~,ind_sort] = sort(b(2,ind_none));
ac(ind_none,:) = ac(ind_none(ind_sort),:);
figure; hold on;
temp = [ac(ind_infra,:); ac(ind_none,:)];
p = pcolor(tac,1:N,temp); set(p,'LineStyle','none')
line([tac(1) tac(end)],[1 1]*Ninfra,'LineWidth',2,'Color','r');
xlim([-900 900]); ylim([1 Ninfra+Nnone]);
caxis([-0.25 0.25]);
title('Automatic Classification');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write the results files (mat)
%close all;
is_infra_all = zeros(1,N); 
is_infra_all(ind_infra) = 1;
mfr_all = mfr;
freq_all = 1./b(2,:);
count0 = 1; count1 = 0;
for n = 1:Nfile
    is_infra = cell(1,Nsheet(n));
    freq = cell(1,Nsheet(n));
    mfr = cell(1,Nsheet(n));
    for m = 1:Nsheet(n)
        count1 = count1+Nunit{n}(m);
        ind = [count0:count1];
        is_infra{m} = is_infra_all([count0:count1]);
        freq{m} = freq_all([count0:count1]);
        mfr{m} = mfr_all([count0:count1]);
        count0 = count1+1;
    end
    save([filepath 'classification\' filename{n} '_res'],'is_infra','freq','mfr');
end


