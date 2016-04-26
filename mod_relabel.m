function [mod2rel] = mod_relabel(mod1,mod2)

% Function matching partition mod2 to a reference partition mod1.
% Currently only works for partitions with the same number of modules,
% and for partitions numbered consecutively from 1 to max(nmod1/2).
%
% Input
% mod1      reference partition
% mod2      partition to be relabeled
%
% Output
% mod2rel   mod2 partition, relabeled for maximum overlap with mod1
%
% Author: Frantisek Vasa (fv247@cam.ac.uk)

% ensure sizes of input match
if all(size(mod1) ~= size(mod2))
    mod1 = mod1';
    if all(size(mod1) ~= size(mod2))
        error('input vectors are of different size')
    end
end

% numbers of modules
nmod1 = max(mod1);
nmod2 = max(mod2);

% ensure that partitions are consecutively numbered
if and(isempty(setdiff(1:nmod1,unique(nmod1))),isempty(setdiff(1:nmod2,unique(nmod2))))
    error('modules in the partitions are not consecutively numbered')
end

% ensure that both partitions have the same number of modules
if nmod1 ~= nmod2
    error('the two partitions have a different number of modules')
end

% evaluate overlap between modules from two partitions
mod_sim = zeros(nmod1,nmod2);
for i = 1:nmod1
    for j = 1:nmod2
        mod_sim(i,j) = sum(and(mod1==i,mod2==j));
    end
end

% re-order overlap matrix so that maximal correlations are lined in
% decreasing order along the diagonal.
x = [];
y = [];

temp_mat = mod_sim;

for i = 1:1:length(mod_sim)
    
    % find max row, col (ignore NaNs, as values which have already been
    % considered get set to NaN.
    [r,c] = find(temp_mat == nanmax(nanmax(temp_mat)));
    x = [x c(1)]; % indexing of first value is used in case of equal overlap
    y = [y r(1)];
    
    temp_mat(r(1),:) = NaN;
    temp_mat(:,c(1)) = NaN;
    
end

% sorted similarity matrix
% mod_sort = mod_sim(y,x);

[b,ix] = sort(x);
y_ord = y(ix);

mod2rel = y_ord(mod2);

end