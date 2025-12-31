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

%%% Polling
presidential.nationalSigmaPoll = 0.0171;
presidential.nationalQPoll = 1.367e-5;

presidential.stateSigmaPoll = 0.0;
presidential.stateQPoll = 0;
presidential.districtSigmaPoll = 0.0304;
presidential.districtQPoll = 1.0635e-5;

%% Generic Ballot Data
%%% Fundamentals Model
genericBallot.nationalMean = 0.5167;
genericBallot.nationalIncumbency = -0.0203;

%%% Uncertainty
genericBallot.nationalSigma = 0.0241;

%%% Polling
genericBallot.nationalSigmaPoll = 0.0167;
% genericBallot.nationalQPoll = 1.4459e-5;
genericBallot.nationalQPoll = 3.9383-6;

%% House Data
%%% Fundamentals Model
house.incumbency = 0.0102;
house.incumbencySigma = 0.0148;
house.presModelSigma = 0.029;
house.prevModelSigma = 0.032;

%%% Polling
house.districtSigmaPoll = 0.0323;
house.districtQPoll = 5.7026e-6;
house.districtIncQPoll = 1.4852e-6;

%% Shared Covariance Data
% Presidential, Generic Ballot
shared.nationalPoll = [1, 0.8876;
    0.8876, 1];

% Shared Covariance
shared.districtSigma = 0.0223;
shared.districtSigmaPoll = 0.0304;
shared.districtQPoll = 3.372e-6;

save("Config.mat")