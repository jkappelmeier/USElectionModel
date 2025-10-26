clear
clc
close all

addpath("..\")
addpath("..\Inputs\")
electionInfo = which("ElectionInfo.csv");
preElectionData = which("PreElectionData.csv");
correlationData = which("CorrelationData.mat");
config = which("Config.mat");

model = Core.Model("2024", electionInfo, preElectionData, correlationData, ...
    config);
xSims = model.simulate(10000);

Visualization.VisualizeResults.printResults(model, xSims);