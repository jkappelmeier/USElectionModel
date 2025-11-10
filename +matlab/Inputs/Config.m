clear

%% Race Specific Strings
electionDate = datetime(2024,11,5); % Election Date
startDate = datetime(2023,11,5); % Campaign Start Date

%% Presidential Data
%%% Fundamentals Model
nationalMean = 0.4926;
incumbency = 0.0238;

%%% Uncertainty
nationalSigma = 0.0449;
stateSigma = 0.0259;
districtSigma = 0.0409;

%%% Polling Uncertainty
nationalSigmaPoll = 0.0171;
nationalQPoll = 1.367e-5;

stateSigmaPoll = 0.0145;
stateQPoll = 9.7683e-6;
districtSigmaPoll = 0.0229;
districtQPoll = 2.4361e-5;

save("Config.mat")