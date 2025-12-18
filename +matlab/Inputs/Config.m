clear

%% Race Specific Strings
electionDate = datetime(2024,11,5); % Election Date
startDate = datetime(2023,11,5); % Campaign Start Date

%% Presidential Data
%%% Fundamentals Model
presidential.nationalMean = 0.4926;
presidential.incumbency = 0.0238;

%%% Uncertainty
presidential.nationalSigma = 0.0449;
presidential.stateSigma = 0.0;
presidential.districtSigma = 0.0396;

%% Generic Ballot Data
%%% Fundamentals Model
genericBallot.nationalMean = 0.5167;
genericBallot.nationalIncumbency = -0.0203;

%%% Uncertainty
genericBallot.nationalSigma = 0.0241;

%% Polling Uncertainty
%%% Presidential
presidential.nationalSigmaPoll = 0.0171;
presidential.nationalQPoll = 1.367e-5;

presidential.stateSigmaPoll = 0.0;
presidential.stateQPoll = 0;
presidential.districtSigmaPoll = 0.0304;
presidential.districtQPoll = 4.3084e-5;

%%% Generic Ballot
genericBallot.nationalSigmaPoll = 0.0167;
genericBallot.nationalQPoll = 1.4459e-5;

%% Cross-Race Correlatio
% Presidential, Generic Ballot
rho.nationalPoll = [1, 0.8876;
    0.8876, 1];

save("Config.mat")