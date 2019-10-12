function [Grps, SingVals] = SpectralClusterings(CKSym,n)
%--------------------------------------------------------------------------
% This function takes an adjacency matrix of a graph and computes the 
% clustering of the nodes using spectral clustering.
% CMat: NxN adjacency matrix
% n: number of groups for clustering
% K: number of largest elements to choose from each column of CMat
% Grps: [grp1,grp2,grp3] for three different types of Spectral Clustering
% SingVals: [SV1,SV2,SV3] singular values for three different types of SC
% LapKernel(:,:,i): last n columns of the Laplacian kernel to apply KMeans
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

%disp('salam')
warning off;
N = size(CKSym,1);
MAXiter = 1000; % Maximum number of iterations for KMeans
REPlic = 20; % Number of replications for KMeans


% Normalized Spectral Clustering according to Ng & Jordan & Weiss
DN = diag( 1./sqrt(sum(CKSym)+eps) );
LapN = speye(N) - DN * CKSym * DN;
[ZZ,sN,vN] = svd(LapN);
kerNS = vN(:,N-n+1:N);
for i = 1:N
    kerNS(i,:) = kerNS(i,:) ./ norm(kerNS(i,:)+eps);
end
%kerNS = DN * kerNS;
SingVals = diag(sN);
Grps = kmeans(kerNS,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');