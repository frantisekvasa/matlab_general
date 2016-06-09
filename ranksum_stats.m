function [r,U] = ranksum_stats(rsum,n1,n2)

% Function to calculate the Mann-Whitney U statistic and the rank-biserial
% correlation based on the sample sizes and the output returned by Matlab
% function "ranksum" - namely, the sum of ranks of the first variable input
% into ranksum.
%
% If running ranksum as [p,h,stats] = ranksum(x,y), then:
% Input:
% rsum      stats.ranksum   (ranksum test statistic)
% n1        length(x)       (size of 1st sample)
% n2        length(y)       (size of 2nd sample)
%
% Output:
% r         rank-biserial correlation()
% U         Mann-Whitney U statistic
%
% Reference:
% Dave S. Kerby (2014) The simple difference formula: an approach to 
% teaching nonparametric correlation. Innovative Teaching 1(3).
% http://www.amsciepub.com/doi/pdf/10.2466/11.IT.3.1

U = rsum-(n1*(n1+1)/2);         % U statistic
r = 1-(2*U)/(n1*n2);            % rank-biserial correlation

end
