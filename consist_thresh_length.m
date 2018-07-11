function [mat_cons,mat_cons_wei] = consist_thresh_length(mat_all,dist,lh_id,rh_id,weighted)

% Function to retain most consistent edges derived from (deterministic)
% tractography, but after estimating consistency within "bins" of edges
% based on their length (to avoid preferential retention of short edges).
% The procedure is separately performed on within- and between-hemispheric
% edges, to avoid preferential retention of one class or the other.
%
% The idea was first proposed in the following paper:
% Misic, Betzel et al. (2015) Cooperative and Competitive Spreading
% Dynamics on the Human Connectome, Neuron 86, 1518–1529.
%
% This was then subsequently applied in:
% Seidlitz et al. (2018) Morphometric Similarity Networks Detect Microscale
% Cortical Organization and Predict Inter-Individual Cognitive Variation, 
% Neuron 97, 231–247.
%
% Frantisek Vasa & Jakob Seidlitz, June 2016 - July 2018 (fv247@cam.ac.uk)
%
% Inputs:
% mat_all           array of tractography-derived (sparse) matrices,
%                   of dimension [N(ROI) x N(ROI) x N(subj)]
%
% dist              distance matrix - here, a single matrix common across
%                   all subjects, derived for example using euclidean
%                   distance between region centroids
%                   (alternatively, the code could be modified to use e.g.
%                   subject-specific distance matrices from tractography)
%
% lh_id             ID's of regions in left hemisphere
% rh_id             ID's or regions in right hemisphere
%                   -> e.g. for a matrix with 100 regions per hemisphere,
%                   ordered as LH -> RH, the ID's would be
%                   lh_id = 1:100; rh_id = 101:200
%
% optional:
% weighted          binary flag determining whether final matrix is weighted
%                   by the median weight (across subjects) of retained edges
%                   (default = 1 == weighted output)
%
% Outputs
% mat_cons         length-based consensus-thresholded matrix (binary)
%
% mat_cons_wei      IF weighted = 1: length-based consensus-thresholded
%                   matrix, weighted by median weight of retained edges

% basic tests
if size(mat_all,1) ~= size(mat_all,2)
    error('incorrect dimensions of connectivity matrix array')
end

if size(dist,1) ~= size(dist,2)
    error('incorrect dimensions of distance matrix')
end

if and(size(mat_all,1) ~= size(dist,1),size(mat_all,2) ~= size(dist,2))
    error('dimensions of connectivity and distance matrices do not match')
end

% default option is to output weighted matrix as well
if nargin < 5
    weighted = 1;
end

nroi = size(mat_all,1);     % number of regions
ns = size(mat_all,3);       % number of subjects

% set diagonal of matrices to zero
for i = 1:1:ns
    for j = 1:1:nroi
        mat_all(j,j,i) = 0;
    end
end

% ID's of within- and between-hemispheric edges
% "hemiw" = within-hemisphere edges
hemiw_id = zeros(nroi); hemiw_id(lh_id,lh_id) = 1; hemiw_id(rh_id,rh_id) = 1;
hemiw_id_all = logical(hemiw_id);                       % full matrix
hemiw_id = logical(and(triu(ones(nroi),1),hemiw_id));   % upper triangular only

% "hemib" = between-hemisphere edges
hemib_id = zeros(nroi); hemib_id(lh_id,rh_id) = 1; hemib_id(rh_id,lh_id) = 1;
hemib_id_all = logical(hemib_id);                       % full matrix
hemib_id = logical(and(triu(ones(nroi),1),hemib_id));   % upper triangular only

% ID matrices can be imaged to verify they look as expected
% figure; imagesc(hemib_id); colormap(flipud(gray));
% figure; imagesc(hemib_id_all); colormap(flipud(gray));
% figure; imagesc(hemiw_id); colormap(flipud(gray));
% figure; imagesc(hemiw_id_all); colormap(flipud(gray));

% calculate edge density within and between hemispheres
for i = 1:1:ns
    temp = mat_all(:,:,i);
    hemiw_dens(i) = sum(temp(hemiw_id)); % within
    hemib_dens(i) = sum(temp(hemib_id)); % between
end

% number of bins within and between hemispheres
% following Misic, Betzel et al. Neuron 2015,
% this is set "heuristically, as the square root of the mean binary density
% across participants", and "this procedure [is] carried out separately
% for inter- and intra-hemispheric edges"
nbins_w = round(sqrt(mean(hemiw_dens(setdiff(1:ns,id_rm))))); % within
nbins_b = round(sqrt(mean(hemib_dens(setdiff(1:ns,id_rm))))); % between

% "distance" matrices (to be thresholded)
dist_w = zeros(nroi); dist_w(hemiw_id) = dist(hemiw_id); % within
dist_b = zeros(nroi); dist_b(hemib_id) = dist(hemib_id); % between

% bin distances into quantile bins (based on number heuristically determined above)
qnt_w = quantile(sort(dist(hemiw_id)),nbins_w-1); qnt_w(end+1) = max(dist(hemiw_id)); % within
qnt_b = quantile(sort(dist(hemib_id)),nbins_b-1); qnt_b(end+1) = max(dist(hemib_id)); % between

% identify "mask" of edges within each bin
% within hemispheres
cut_l = 0; % lower cut-off - intially set to 0, then iteratively re-defined as code proceeds into "longer" bins
dist_w_bin = zeros(nroi,nroi,nbins_w);
for b = 1:nbins_w
    
    temp = zeros(nroi);
    temp(and(dist_w>cut_l,dist_w<qnt_w(b))) = 1; % find edges whose distance fits within limits of the current bin
    dist_w_bin(:,:,b) = temp; % store "mask" of those edges
    
    cut_l = qnt_w(b);
end

% between hemispheres
cut_l = 0; % lower cut-off - intially set to 0, then iteratively re-defined as code proceeds into "longer" bins
dist_b_bin = zeros(nroi,nroi,nbins_b);
for b = 1:nbins_b
    
    temp = zeros(nroi);
    temp(and(dist_b>cut_l,dist_b<qnt_b(b))) = 1; % find edges whose distance fits within limits of the current bin
    dist_b_bin(:,:,b) = temp; % store "mask" of those edges
    
    cut_l = qnt_b(b);
end

cons_mat_b = sum(logical(mat_all,3)); % matrix of (binary) edge consistency

% perform binning
% as in Misic, Betzel et al Neuron 2015: "if the mean number of edges
% (across participants) in a particular bin i is equal to k_i, we selected
% the k_i most commonly occurring edges in that bin"
%
% within bin
edges_w_bin = zeros(nroi,nroi,nbins_w);
for b = 1:1:nbins_w
    
    temp = cons_mat_b.*dist_w_bin(:,:,b);       % consistency of edges in current bin
    nedge_bin = round(sum(temp(:))/ns);         % mean number of edges in current bin
    temp_sort = sort(temp(temp>0),'descend');   % sort edges by consistency
    temp_mat_bin = zeros(nroi);
    temp_mat_bin(temp>=temp_sort(nedge_bin)) = 1;   % retain most consistent edges (equal in number to mean across participants)
    edges_w_bin(:,:,b) = temp_mat_bin;              % store edges
    
end
mat_w = sum(edges_w_bin,3); % "collapse" within-hemi matrix across bins

% between bin
edges_b_bin = zeros(nroi,nroi,nbins_b);
for b = 1:1:nbins_b
    
    temp = cons_mat_b.*dist_b_bin(:,:,b);       % consistency of edges in current bin
    nedge_bin = round(sum(temp(:))/ns);         % mean number of edges in current bin
    temp_sort = sort(temp(temp>0),'descend');   % sort edges by consistency
    temp_mat_bin = zeros(nroi);
    temp_mat_bin(temp>=temp_sort(nedge_bin)) = 1;   % retain most consistent edges (equal in number to mean across participants)
    edges_b_bin(:,:,b) = temp_mat_bin;              % store edges
    
end
mat_b = sum(edges_b_bin,3); % "collapse" between-hemi matrix across bins

% the final matrix is a combination of within- and between-hemisphere bins
% as well as the sum of itself and its transpose (as the above operations
% were performed on the upper triangular only)
mat_cons = mat_w+mat_b;
mat_cons = mat_cons+mat_cons';

% optional weighted version
if weighted
    % median of non-zero entries for weighted matrix
    temp_mat_all_wei = mat_all; temp_mat_all_wei(temp_mat_all_wei==0) = NaN; % first, set zero entries to NaN
    median_mat_wei = nanmedian(temp_mat_all_wei,3); % then apply nanmedian to ignore these (initially zero) entries
    % weight retained edges by median weights
    mat_cons_wei = zeros(nroi);
    mat_cons_wei(logical(mat_cons)) = median_mat_wei(logical(mat_cons));
end

end
