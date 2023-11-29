function [covarmat, precmat] = stable_partcorr(adjmat, samplecov)

%MATLAB implementation of the Maximum Likelihood Estimate of a covariance
%matrix, given the sample covariance matrix and a fixed Gaussian graphical
%model topology, as outlined in:
% The Elements of Statistical Learning: Data Mining, Inference, and
% Prediction by T. Hastie, R. Tibshirani, and J. Friedman, 
% Chapter 17.3.1 "Estimation of the Parameters when the Graph Structure is Known"
% available at https://doi.org/10.1007/978-0-387-84858-7_17 

%adjmat is a logical square matrix indicating the adjacency matrix of the
%undirected graphical model topology.
%samplecov is the sample covariance matrix of the data being modeled.
%adjmat & samplecov must both be n x n matrices for the same n


    n = size(adjmat,1);
    adjmat(logical(eye(n))) = 0;
    
    precmat = zeros(size(samplecov));
    W = samplecov;
    W_iter = samplecov; %initialize
    tol = 1e-6; %convergence threshold
    converg_flag = 0;
    count = 0;
    
    while converg_flag == 0
        count = count+1;
        for j = 1:n
            idx = true(1,n); idx(j) = false;
            W_11 = W_iter(idx,idx);
            beta = zeros(n-1,1);
            neighb_id = adjmat(idx,j);
            
            W_11s = W_iter(idx & adjmat(j,:), idx & adjmat(j,:));
            s_12s = samplecov(idx & adjmat(j,:),j);
            betas = W_11s\s_12s;
            
            beta(neighb_id) = betas;
            
            W_iter(j,idx) = (W_11*beta)';
            W_iter(idx,j) = W_11*beta;
        end
        
        if max(abs(W(:)-W_iter(:)))<tol
            covarmat = W_iter;
            precmat = inv(covarmat);
            converg_flag=1;
        else
            W = W_iter;
        end
    end
    disp(count)
end