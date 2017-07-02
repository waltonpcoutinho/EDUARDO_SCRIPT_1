function [sicd idx fconv] = kmeans(data_path,qtd,dim,K);

% k-means spike sorting 
% 
% K: number of classes
% X: dataset matrix (num. of cases x dimensions)
%
% cluster = k_means(K,X);
%
% Written by Sturla Molden, sturla@molden.net
	tic;

	range=[0 0 qtd-1 dim-1];
    X = dlmread(data_path,',',range);

	[n,d] = size(X);
	% Start with random assignments 
	clust = ones(n,1);
	for t = 1:n;
	  	randomsample = randperm(K);
	  	clust(t) = randomsample(1);
	end

	% Set parameters
	clust2 = clust;
	g = max(clust);
	term = 1;

	best_sicd = inf;
	nImp = 1;

	% Run K-means until convergence
	%disp('Sorting spikes ...');
	while term ~= 0;

		% Find centroids
		centroids = [];
		for c = 1:g;
			index = find(clust == c);
			if isempty(index) ~= 1; 
				if length(index) == 1
					centroids = [centroids; X(index,:)];
				else
					centroids = [centroids; mean(X(index,:))];
				end
			end
		end 

		% Find distances to centroids and recluster 
		[g,dim] = size(centroids);
		dist = ones(n,1);
		for s = 1:n;
		 	x = X(s,:); 
			x = ones(g,1)*x;
			d = (centroids - x).^2;
			d = sqrt(sum(d')');
			[m,index] = min(d);	
			clust2(s) = index;
			dist(s) = m; 
		end 

		% Check for convergence
		term = sum(clust ~= clust2); 
		clust = clust2; 

		sicd = sum(dist);

		%f(reshape(centroids',1,[]), K*dim, K, X, qtd)

		if sicd < best_sicd
			idx = clust2;
			cents = centroids';
			cent_best = cents(:);
			cent_best = cent_best';
			best_sicd = sicd;
			fconv(nImp,1) = toc;
			fconv(nImp,2) = sicd;
			nImp = nImp + 1;
		end

	end

end

