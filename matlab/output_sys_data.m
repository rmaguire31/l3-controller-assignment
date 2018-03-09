function output_sys_data(sys, dirout, suffix)

mkdir(dirout);

% Write root locus to file for use in pgf plots
[r,k] = rlocus(sys);
rloc = ones(size(k,2), size(k,1) + 2*size(r,1));
rloc(:,1:size(k,1)) = k';
header = 'K';
for j = 1:size(r,1)
    header = sprintf('%s\tre%d', header, j);
    rloc(:,size(k,1)+2*j-1) = real(r(j,:))';
    header = sprintf('%s\tim%d', header, j);
    rloc(:,size(k,1)+2*j) = imag(r(j,:))';
end
fname = [dirout '/rlocus' suffix '.dat'];
fid = fopen(fname, 'wt');
fprintf(fid, '%s\n', header);
fclose(fid);
dlmwrite(fname, rloc, '-append', 'delimiter', sprintf('\t'));

% Write step responses to file for use in pgfplots
t = 0:1e-2:20;
stepres = ones(length(t),2);
stepres(:,1) = t';
header = sprintf('t\ty');
y = step(feedback(sys, 1), t);
stepres(:,2) = y';
fname = [dirout '/stepres' suffix '.dat'];
fid = fopen(fname, 'wt');
fprintf(fid, '%s\n', header);
fclose(fid);
dlmwrite(fname, stepres, '-append', 'delimiter', sprintf('\t'));

stepinfo(feedback(sys, 1))