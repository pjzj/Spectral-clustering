function [NI,sigma,dist2]=Knearest_4(data,k,ind)
data=single(data');
[N,D]=size(data');
NI = cell(1,N);
if ind==2
    if D>N
        a = sum(data.*data);
        dist2 = sqrt(repmat(a',[1 N]) + repmat(a,[N 1]) - 2*(data'*data));
        for i=1:N
            [dist_sort,J]=sort(dist2(:,i));
            Ii=J(1:k+1);
            NI{i} = Ii;
        end;
        sigma=mean(mean(dist2));
    else
        temp=[];
        for i=1:N
            x = data(:,i);
            dist2 = sum((data-repmat(x,[1 N])).^2,1);
            [dist_sort,J] = sort(dist2);
            Ii = J(1:k+1);
            NI{i} = Ii;
            temp=[temp;dist2];
        end;
        dist2=temp;
        clear temp;
        sigma=mean(mean(dist2));
    end
end;

if ind==1
    temp=[];
    for i=1:N
        x = data(:,i);
        dist2 = sum(abs(data-repmat(x,[1 N])),1);
       [dist_sort,J] = sort(dist2);
       Ii = J(1:k+1);
       NI{i} = Ii;
        temp=[temp;dist2];
    end;
    dist2=temp;
    clear temp;
    sigma=mean(mean(dist2));
end
% if ind==1
%     temp=[];
%     for i=1:N
%         x = data(:,i);
%         dist2 = sum(abs(data-repmat(x,[1 N])),1);
%         [dist_sort,J] = sort(dist2);
%         Ii = J(1:k+1);
%         NI{i} = Ii;
%         temp=[temp;dist2];
%     end;
%     dist2=temp;
%     clear temp;
%     sigma=mean(mean(dist2));
% end
if ind==3
    temp=[];
    for i=1:N
        data(:,i) = data(:,i)/norm(data(:,i));
    end
    for i=1:N
         x = data(:,i);
        dist2 = sum(data.*repmat(x,[1 N]),1);
        [dist_sort,J] = sort(dist2,'descend');
        Ii = J(1:k+1);
        NI{i} = Ii;
        temp=[temp;dist2];
    end;
    dist2=temp;
    clear temp;
    sigma=mean(mean(dist2));
end
if ind==4
    temp=[];
   for i=1:N
            x = data(:,i);
            dist22 = sum((data-repmat(x,[1 N])).^2,1);
            dist21= sum(abs(data-repmat(x,[1 N])),1);
            dist2=dist21+dist22;
            temp=[temp;dist2];
   end

   parfor i=1:N
        data(:,i) = data(:,i)/norm(data(:,i));
    end
    parfor i=1:N
            dist23 = sum(data.*repmat(x,[1 N]),1);
            temp(i,:)=temp(i,:)+acos(dist23);
            [dist_sort,J] = sort(temp(i,:));
            Ii = J(1:k+1);
            NI{i} = Ii;            
    end
        dist2=temp;
       clear temp; 
       sigma=mean(mean(dist2));
end

if ind==5
    temp=[];
    for i=1:N
        x = data(:,i);
        dist22 = sum((data-repmat(x,[1 N])).^2,1);
        dist21= sum(abs(data-repmat(x,[1 N])),1);
        dist2=dist21+dist22;
        [dist_sort,J] = sort(dist2);
        Ii = J(1:k+1);
        NI{i} = Ii;
        temp=[temp;dist2];
    end
    dist2=temp;
    clear temp;
    sigma=mean(mean(dist2));
end 
