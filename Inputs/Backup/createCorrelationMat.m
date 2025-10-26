clear

corrMatrix = readmatrix("CorrelationMatrix.csv");
corrInfo = readtable("CorrelationInfoData.csv", "TextType", "string");

save("CorrelationData.mat");

