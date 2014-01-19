function [ acc, pred ] = svm_ovr_class( trainData, trainLabel, testData, testLabel )
%SVM_OVR_PRED SVM one-vs-rest classification.

numLabels=max(trainLabel);
model = cell(numLabels,1);
for i=1:numLabels
    model{i} = svmtrain(double(trainLabel==i), trainData, '-q -c 100 -t 0 -b 1');
end

% get probability estimates of test instances using each model
prob = zeros(size(testData,1),numLabels);
for i=1:numLabels
    [~,~,p] = svmpredict(double(testLabel==i), testData, model{i}, '-b 1');
    prob(:,i) = p(:,model{i}.Label==1);    %# probability of class==k
end

% predict the class with the highest probability
[~,pred] = max(prob,[],2);
acc = sum(pred == testLabel) ./ numel(testLabel);    %# accuracy
%C = confusionmat(testLabel, pred)                   %# confusion matrix

end

