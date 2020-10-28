function[x,text,sheet] = load_data_review(filename)
filepath = 'Data\';
%fileinfo
[~,sheet] = xlsfinfo([filepath filename]);
Nsheet = numel(sheet);
%load neurons
count = 0;
x = []; text = cell(Nsheet,100);
for n = 1:Nsheet
    [x_temp, text_temp] = xlsread([filepath filename],sheet{n});
    x{n} = x_temp(:,1:end);
    text(n,1:size(text_temp,2)) = text_temp(1,:);
    disp(sprintf('loaded %s expts out of %s',num2str(n),num2str(Nsheet)));
end
save([filepath filename '_data'],'x','text','sheet','Nsheet');