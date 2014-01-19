fprintf('AP:\t');
for i = 1:length(results.res{1})
    fprintf('%f\t',results.res{1}(i).AP*100);
end
fprintf('\n');

fprintf('prec@10:\t');
for i = 1:length(results.res{1})
    fprintf('%f\t',results.res{1}(i).prec10*100);
end
fprintf('\n');

fprintf('prec@50:\t');
for i = 1:length(results.res{1})
    fprintf('%f\t',results.res{1}(i).prec50*100);
end
fprintf('\n');