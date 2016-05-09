function [A MST] = mst_threshold(W,dens,bin)

% Function to threshold weighted matrices, enforcing connectedness by 
% setting the MAXIMUM spanning tree as a backbone to the network.
% Works only with symmetric matrices. Ignores negative weights.
%
% Dependence: matlab_bgl (at http://dgleich.github.io/matlab-bgl/)
% Specifically, kruskal_mst (which itself depends on other matlab_bgl
% functions).
%
% Input
% W     weighted matrix
% dens  density (in range [0,1])
% bin   "binary" = optional flag for output: 
%       T if binary (default)
%       F if weighted
%
% Output
% A     thresholded matrix
% MST   the maximum spanning tree
%
% author: Petra Vertes
% modified by Frantisek Vasa (fv247@cam.ac.uk)

if nargin < 3
    bin = true;
end

nnod = size(W,1);           % N nodes
W = (W+W')/2;               % force symmetrize matrix
W(W<0) = 0;                 % set negative Wrrelations to 0
W(1:nnod+1:nnod*nnod) = 1;  % set diagonal to ones

% create MST (the MAXIMUM spanning tree of the network)  
% kruskal_mst finds the minimum spanning tree, but it is here applied to
% the inverse of the weighted input matrix, W.
D=ones(nnod,nnod)./W;
MST=kruskal_mst(sparse(D));

% order C according to decreasing weights in the Wrrelation matrix
W_tr = triu(W,1);
ind = find(triu(ones(nnod),1));
Wlist = W_tr(ind);
[~,IX] = sort(Wlist,'descend');
[row,col]=ind2sub([nnod,nnod],ind(IX));

% store initial MST in the adjacency matrix A that defines the output network
if bin
    A = double(logical(full(MST)));
else
    A = zeros(nnod);
    A(logical(full(MST))) = W(logical(full(MST)));
end

MST = A; % store final MST (for 2nd output)

% grow the network according to weights in W matrix 
t = 1;
enum = nnod-1;
% add edges in correct order until all possible edges exist
while (enum < dens*nnod*(nnod-1)/2)
    % if edge wasn't initially included in MST
    if A(row(t),col(t)) == 0
        if bin % binary version
            A(row(t),col(t)) = 1; 
            A(col(t),row(t)) = 1;
            enum=enum+1;
        else % weighted version
            A(row(t),col(t)) = W(row(t),col(t)); 
            A(col(t),row(t)) = W(col(t),row(t));
            enum=enum+1;
        end
    end
    t=t+1;
end

end
