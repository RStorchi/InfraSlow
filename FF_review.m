function[] = FF_review()
%
count = 1;
filename = [];
filename{count} = 'dLGN_RC'; count = count+1;
filename{count} = 'dLGN_MELKO'; count = count+1;
filename{count} = 'dLGN_RDCL'; count = count+1;

filename{count} = 'vLGN_RC'; count = count+1;
filename{count} = 'vLGN_MELKO'; count = count+1;
filename{count} = 'vLGN_RDCL'; count = count+1;

filename{count} = 'OPN_RC'; count = count+1;
filename{count} = 'OPN_MELKO'; count = count+1;
filename{count} = 'OPN_RDCL'; count = count+1;

filename{count} = 'pret_RC'; count = count+1;
filename{count} = 'pret_MELKO'; count = count+1;
filename{count} = 'pret_RDCL'; count = count+1;
%
Nfile = numel(filename); 
filepath1 = 'Data\';
filepath2 = 'Data\classification\';

%
infra = []; xb = []; 
for n = 1:Nfile
    load([filepath1 filename{n} '_data'],'x');
    load([filepath2 filename{n} '_infra_calc_res'],'is_infra');
    xb = [xb x];
    infra = [infra horzcat(is_infra{:})];
end
N = numel(infra);
%
xmax = 900; count = 0;
for n = 1:numel(xb)
    for m = 1:size(xb{n},2)
        count = count+1;
        temp = xb{n}(:,m)';
        temp = temp(temp>0);
        [FF(count,:),T] = FF_review_v2_sub1(temp,xmax);
        disp(sprintf('unit %s of %s',num2str(count),num2str(N)));
    end
end
%
save([filepath1 'FF_calc'],'T','FF','infra');





function[FF,T] = FF_review_v2_sub1(x,xmax)
T = 10.^[-3:0.25:2];
NT = numel(T);
for n = 1:NT
    Tvec = [0:T(n):xmax];
    xbin = hist(x,Tvec);
    FF(n) = var(xbin)/mean(xbin);
end