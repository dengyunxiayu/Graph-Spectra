clc; clear ; close all;
E = csvread('example1.dat');
col1 = E(:,1);
col2 = E(:,2);
max_ids = max(max(col1,col2));
G = graph(col1,col2);

figure('Name','Original graph without clusters');
p = plot(G);
p.LineWidth= 0.2500;
p.LineStyle= '-';
p.Marker = '.';
p.EdgeColor = 'k';

As= sparse(col1, col2, 1, max_ids, max_ids); 
A = full(As);
Di = diag(sum(A,2));
L = (Di^(-0.5)) * A * (Di^(-0.5));

[X,D1] = eigs(L,size(L,1));
eigenvalues = diag(D1);
eig_gaps = -1*diff(eigenvalues); 
array_k = diff(sort(eigenvalues,'descend'));
[w,k]=min(array_k); 
X = X(:,(1:k));
Y = X./(sum(X.^2,2)).^(0.5);
Idx = kmeans(Y, k);
figure('Name','Graph with clusters');
p = plot(G);
p.EdgeColor = 'k';
p.LineWidth= 0.2500;
p.LineStyle= '-';
p.Marker = '.';
p.MarkerSize = 10;
cmap = lines();

for i=1:size(Idx,1)
    cluster = Idx(i,1);
    highlight(p,i,'NodeColor',cmap(cluster,:));
end

figure('Name', 'Sparsity Pattern');
spy(A);

new_L=sparse(Di-A);
[X2,new_DD] = eigs(new_L,k,'smallestreal');
figure('Name','Sorted Fiedler Vector');
plot(sort(X2(:,2)));