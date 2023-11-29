function [stability_matrix, LLH_matrix] = stable_gLASSO(DataFull,varid,lambda_list,k,n)
    
    % this code depends on the gLASSO MATLAB implementation that can be
    % found at: https://github.com/xiaohuichen88/Graphical-Lasso 
    addpath 'matlab_mtl'
%%
    %DataFull is MATLAB table object containing the data to be modeled
    
    %varid is an index vector indicating the columns of DataFull to use for
    %analysis

    %lambda_list is an ordered list containing the LASSO lambda parameters
    %on which to evaluate graphical LASSO

    % k is k-fold cross-validation parameter (i.e. what proportion of data
    % to leave out during subsampling

    % n is number of data subsamples to take
%%

    data = preprocess_data(DataFull{:,varid}); %data-matrix to be analyze for gLASSO
    
    p = size(data,2);
    lambda_res = numel(lambda_list);
    zero_cutoff = 1e-4;
    
    LLH_matrix = zeros(n,lambda_res); %predictive power at different lambdas, from cross-valid
    stability_matrix = zeros(p^2,lambda_res); %percentage of runs at given lambda that include a given edge as non-zero
    
    for i = 1:lambda_res
        stability_counter = zeros(p^2, n);
        parfor j = 1:n
            [testdata,traindata] = splitdata(data,k); %split data k-fold into training and test
            [Theta_train,W_train] = graphicalLasso(cov(traindata),lambda_list(i));%run gLASSO on training
            stability_counter(:,j) = abs(Theta_train(:))>zero_cutoff;
            LLH_matrix(j,i) = -gaussianLLH(testdata,W_train); %evaluate on test; save results LLH_matrix  
        end
        stability_matrix(:,i) = mean(stability_counter,2);%evaluate statistics for stability_counter
    end
end

%% auxilary functions

function data = preprocess_data(data_in)
    null_id = any(isnan(data_in),2);
    data = zscore(data_in(~null_id,:));
end

function [testdata, traindata] = splitdata(data,k)
    n_sample = size(data,1);
    test_id = randperm(n_sample,floor(n_sample/k));
    train_id = setdiff(1:n_sample, test_id);
    testdata = data(test_id,:);
    traindata = data(train_id,:);
end

function LLH = gaussianLLH(sampledata, covarmat)
    %sampledata must come from a distribution with mean-zero
    LLH = sum(log(mvnpdf(sampledata,zeros(numel(sampledata,2),1),covarmat)));
end