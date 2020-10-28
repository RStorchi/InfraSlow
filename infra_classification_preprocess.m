function[] = infra_classification_preprocess(filename)
%load data
filepath = 'Data\';
load([filepath filename '_data']);
%bin data
dT = 1;
Tbin = 0:dT:900; Nt = numel(Tbin);
xbin = cell(1,Nsheet);
xac =  cell(1,Nsheet);
xmfr = cell(1,Nsheet);
xok = cell(1,Nsheet);
N = 0;
for n = 1:Nsheet
    for m = 1:size(x{n},2)
        xbin{n}(:,m) = hist(x{n}(:,m),Tbin);     
        [xac{n}(:,m),tac] = xcorr(detrend(xbin{n}(:,m)),'coeff');
        xmfr{n}(m) = mean(xbin{n}(:,m))/dT;
        xok{n}(m) = sum(isnan(xac{n}(:,m)))==0;
        N = N+1;
    end
end
tac = tac';
ind_pos = find(tac>0);
%fit autocorr
errall0 = cell(1,Nsheet);
ball1 = cell(1,Nsheet);
errall1 = cell(1,Nsheet);
ball2 = cell(1,Nsheet);
errall2 = cell(1,Nsheet);
xac_fit = cell(1,Nsheet);
is_infra = cell(1,Nsheet);
ball3 = cell(1,Nsheet);
%%%%%%%%%%%%%%%%%%%%%%%infra criterion%%%%%%%%%%%%%%%%%%%%%%%
THper = [50 450]; %period
THdec = 1; %decay
THerr = 0.15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 0;
for n = 1:Nsheet
    for m = 1:size(x{n},2)
        if xok{n}(m)
            errall0{n}(m) = norm(xac{n}(ind_pos,m)-mean(xac{n}(ind_pos,m)));
            b0 = [0.1 100]; 
            b = fmincon(@(b) xcorr_fun1(b,tac(ind_pos),xac{n}(ind_pos,m)),b0,[],[],[],[],[0 1],[1 18000]);
            errall1{n}(m) = xcorr_fun1(b,tac(ind_pos),xac{n}(ind_pos,m)); 
            ball1{n}(:,m) = b;
            b0 = [0.1 100 100 0.1 200 200]; 
            b = fmincon(@(b) xcorr_fun2(b,tac(ind_pos),xac{n}(ind_pos,m)),b0,[],[],[],[],[0 1 1 0 1 1],[1 900 1800 1 900 18000]);
            errall2{n}(m) = xcorr_fun2(b,tac(ind_pos),xac{n}(ind_pos,m)); 
            ball2{n}(:,m) = b;        
            %
            xac_fit{n}(:,m) = b(1)*cos(2*pi*tac(ind_pos)/b(2)).*exp(-tac(ind_pos)/b(3))+b(4)*cos(2*pi*tac(ind_pos)/b(5)).*exp(-tac(ind_pos)/b(6));
        else
            ball1{n}(:,m) = [NaN; NaN];
            ball2{n}(:,m) = [NaN; NaN; NaN; NaN; NaN; NaN];
            errall1{n}(:,m) = NaN;
            xac_fit{n}(:,m) = NaN*ones(numel(ind_pos),1);
        end
        %  
        count = count+1;
        disp(sprintf('cell %s out of %s',num2str(count),num2str(N)));
        %
        is_graph = false;
        if is_graph %& is_infra{n}(m)
            fig = figure; 
            set(fig,'Position',[50 50 1000 600]);
            subplot(2,2,1), hold on; 
            plot(tac(ind_pos),xac{n}(ind_pos,m),'k',tac(ind_pos),xac_fit{n}(:,m),'b','LineWidth',2);
            title([errall0{n}(m)' errall1{n}(m)' errall2{n}(m)' 100*derr']);
            subplot(2,2,2), hold on; 
            plot(filter(ones(1,3)/3,1,frtemp),'k','LineWidth',2);
            plot(Tbin,c(1)*cos((2*pi*Tbin)/bdom(2)+c(2)),'b','LineWidth',2);
            title(errc);
            subplot(2,2,3), hold on; 
            plot(tac(ind_pos),b(1)*cos(2*pi*tac(ind_pos)/b(2)).*exp(-tac(ind_pos)/b(3)),'r','LineWidth',2);
            plot(tac(ind_pos),b(4)*cos(2*pi*tac(ind_pos)/b(5)).*exp(-tac(ind_pos)/b(6)),'g','LineWidth',2);
            title([b(3)/b(2) b(6)/b(5) is_infra{n}(m)]);
            subplot(2,2,4), hold on; 
            bar(1,b(1)*b(3),'BarWidth',0.5,'FaceColor','r')
            bar(2,b(4)*b(6),'BarWidth',0.5,'FaceColor','g')
            xlim([-2 5]);
            ginput(); close all;
        end
    end
end
%save
save([filepath filename '_infra_calc'],'xbin','xac','xac_fit','tac','is_infra','ball1','errall1','ball2','errall2','THper','THdec','xmfr');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[err] = xcorr_fun1(b,x,y)
yfit = b(1)*exp(-x/b(2));
err = norm(y-yfit);

function[err] = xcorr_fun2(b,x,y)
yfit = b(1)*cos(2*pi*x/b(2)).*exp(-x/b(3))+b(4)*cos(2*pi*x/b(5)).*exp(-x/b(6));
err = norm(y-yfit);


