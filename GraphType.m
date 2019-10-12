function [L,index,X]=GraphType(D,N,Graph_type,k,NI,dist2)
if ~exist('NI', 'var')
    NI = {};
end
if ~exist('sigma', 'var')
    sigma = 1;
end
if ~exist('dist2', 'var')
    dist2 = [];
end
%%
% labels=data(:,end);
% data=data(:,1:end-1);
%¿≠∆’¿≠ÀπÕº
if strcmp(Graph_type,'Laplace')
    %      data=data/max(max(abs(data)));
    %      W=constructW(data,options);
    W=zeros(N,N);
    for i=1:N
        W(i,NI{i}(2:k+1))=exp(-dist2(i,NI{i}(2:k+1))/2);
    end
    W=max(W,W');
%         figure,
%     spy(W)
    D_inv = diag(full(sum(W,2)).^(-0.5));
    L=D_inv*W*D_inv;
    % save laplace L
    %     figure,
    %     spy(L)
end
%%
%LLEÕº
if strcmp(Graph_type,'LLE')
    data=data';
    tol = 1e-4;
    W = zeros(N, N);
    for i=1:N
        Ii = NI{i};%ni(i,2:end);
        Ii=Ii(2:end);
        z = data(:,Ii) - repmat(data(:,i), 1, k);
        C = z' * z;
        C = C + eye(k, k) * tol * trace(C);
        W(Ii,i) = C \ ones(k, 1);
        W(Ii,i) = W(Ii,i) / sum(W(Ii,i));
    end
    W=max(W,W');
    D_inv = diag(full(sum(W,2)).^(-0.5));
    L=D_inv*W*D_inv;
end



%%
if strcmp(Graph_type,'Euclidean')
    X=data;
    n = size(data, 1);
    %[D, ni] = find_nn(X, k+1);
    W=zeros(n,n);
    for i=1:n
        Ii = NI{i};%ni(i,2:end);
        Ii=Ii(2:end);
        for j=1:k
            W(i,Ii(j))=1/norm(X(Ii(j),:));
        end
    end
    W=max(W,W');
    D_inv = diag(full(sum(W,2)).^(-0.5));
    L=D_inv*W*D_inv;
end





