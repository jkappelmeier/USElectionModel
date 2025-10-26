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

save("Config.mat")