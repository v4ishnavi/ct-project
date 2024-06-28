function threshold= iterative_threshold(X, nc)
%here threshold is the output array 
%X is the distribution 
%nc -> number of clusters: here 4 
if isrow(X)
    X = X';
end
threshold = zeros(1, nc);
[~, C] = kmeans(X, nc);
% disp(C)
for p = 1:nc
threshold(p) = C(p,1);
end

end