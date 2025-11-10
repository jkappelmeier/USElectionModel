clear
clc
close all

% Files
electionInfo = "ElectionInfo.csv";
preElectionData = "PreElectionData.csv";
correlationData = "CorrelationData.mat";
config = "Config.mat";
pollsFile = "Polls.csv";

% Setup Fundamentals
model = matlab.Core.Model("2024", electionInfo, preElectionData, ...
    correlationData, config);

% Add Polls
model.addPolls(pollsFile);
model.runPollingAverage();

% Combine Polls and Fundamentals
model.generateEstimate();

% Simulate Elections
xSims = model.simulate(10000);

% Display Results
matlab.Visualization.VisualizeResults.printResults(model, xSims);