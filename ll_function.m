function [loglikelihood] = ll_function(data1,theta)
[testdata] = datagen(data1(:,1),theta);
[testdata_2] = datagen_virus(data1(:,1),theta);
loglikelihood_1  = sum(log(pdf('Normal',log(data1(:,2)),log(testdata'),1)));

loglikelihood_2  = sum(log(pdf('Normal',log(data1(:,3)),log(testdata_2'),1)));

loglikelihood = loglikelihood_1+loglikelihood_2;

% temple_data = log(data(:,2))-log(testdata');
% loglikelihood  = -sum(temple_data.^2);
end