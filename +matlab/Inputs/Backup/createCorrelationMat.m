clear

corrMatrix = readmatrix("CorrelationMatrix.csv");
districtInfo = readtable("CorrelationInfoData.csv", "TextType", "string");
redistrictingData = readtable("RedistrictingData.csv", "TextType", "string");

save("DistrictData.mat");

