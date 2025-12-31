clear

corrMatrix = readmatrix("CorrelationMatrix.csv");
districtInfo = readtable("CorrelationInfoData.csv", "TextType", "string");

save("DistrictData.mat");

