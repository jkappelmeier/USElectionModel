clear

%% Race Specific Strings
electionDate = datetime(2024,11,5); % Election Date
startDate = datetime(2023,11,5); % Campaign Start Date

%% Shared Covariance Data
% Presidential, Generic Ballot, House, Senate, Gubernatorial
shared.nationalFund = [1, 0, 0.4912, 0.5324, 0;
    0, 1, 0.8019, 0.6388, 0.8365;
    0.4912, 0.8019, 1, 0.8769, 0.4718;
    0.5324, 0.6388, 0.8769, 1, 0.2894;
    0, 0.8365, 0.4718, 0.2894, 1];
shared.nationalPoll = ...
    [1, 0.8876, 0.8876, 0.9909, 0.9104;
    0.8876, 1, 1, 0.8459, 0.7934;
    0.8876, 1, 1, 0.8459, 0.7934;
    0.9909, 0.8459, 0.8459, 1, 0.8979;
    0.9104, 0.7934, 0.7934, 0.8979, 1];

% House, Senate, Governor
shared.incCorr = [1, 0.2281, 0.1799;
    0.2281, 1, 0.3960;
    0.1799, 0.3960, 1];

% Shared Covariance
shared.districtSigma = 0.0261;
shared.districtSigmaPoll = 0.025;
shared.districtQPoll = 3.372e-6;

%% Presidential Data
%%% Fundamentals Model
presidential.nationalMean = 0.4926;
presidential.incumbency = 0.0238;

%%% Uncertainty
presidential.nationalSigma = 0.0449;
presidential.stateSigma = 0.0;
presidential.districtSigma = sqrt(0.0396^2 - shared.districtSigma^2);

%%% Polling
presidential.nationalSigmaPoll = 0.0171;
presidential.nationalQPoll = 1.367e-5;

presidential.stateSigmaPoll = 0.0;
presidential.stateQPoll = 0;
presidential.districtSigmaPoll = sqrt(0.0259^2 - shared.districtSigmaPoll^2);
presidential.districtQPoll = presidential.districtSigma^2/presidential.nationalSigma^2*presidential.nationalQPoll - shared.districtQPoll;

%% Generic Ballot Data
%%% Fundamentals Model
genericBallot.nationalMean = 0.5167;
genericBallot.nationalIncumbency = -0.0203;

%%% Uncertainty
genericBallot.nationalSigma = 0.0241;

%%% Polling
genericBallot.nationalSigmaPoll = 0.0167;
% genericBallot.nationalQPoll = 1.4459e-5;
genericBallot.nationalQPoll = genericBallot.nationalSigma^2/presidential.nationalSigma^2*presidential.nationalQPoll;

%% House Data
%%% Fundamentals Model
house.nationalBiasIncumbency = -0.0147;
house.nationalBiasMean = 0.0124;
house.nationalBiasSigma = 0.0207;

house.incumbency = 0.012;
house.incumbencySigma = 0.0148;
house.presModelSigma = sqrt(0.029^2 - shared.districtSigma^2);
house.prevModelSigma = sqrt(0.032^2 - shared.districtSigma^2);

%%% Polling
house.nationalSigmaPoll = genericBallot.nationalSigmaPoll;
house.nationalQPoll = genericBallot.nationalQPoll;

house.districtSigmaPoll = sqrt(0.0307^2 - shared.districtSigmaPoll^2);
house.districtQPoll = house.presModelSigma^2/presidential.nationalSigma^2*presidential.nationalQPoll - shared.districtQPoll;
house.districtIncQPoll = house.incumbencySigma^2/presidential.nationalSigma^2*presidential.nationalQPoll;

%% Senate Data
%%% Initial Information
senate.seats.dem = 28;
senate.seats.rep = 38;

%%% Fundamentals Model
senate.nationalBiasIncumbency = -0.0188;
senate.nationalBiasMean = 0.0222;
senate.nationalBiasSigma = 0.0207;

senate.incumbency = 0.0302;
senate.presPrevSigma = 0.0218;
senate.presModelSigma = sqrt(0.0535^2 - 0.6043 * shared.districtSigma^2 - senate.presPrevSigma^2);
senate.prevModelSigma = sqrt(0.0607^2 - 0.6043 * shared.districtSigma^2 - senate.presPrevSigma^2);
senate.incumbencySigma = 0.014;

%%% Polling
senate.nationalSigmaPoll = 0.016;
senate.nationalQPoll = senate.nationalBiasSigma^2/presidential.nationalSigma^2*presidential.nationalQPoll;

senate.stateSigmaPoll = sqrt(0.0262^2 - 0.6043 * shared.districtSigmaPoll^2);
senate.stateQPoll = senate.presModelSigma^2 / presidential.nationalSigma^2 * presidential.nationalQPoll;
senate.stateIncQPoll = senate.incumbencySigma^2 / presidential.nationalSigma^2 * presidential.nationalQPoll;

%% Gubernatorial Data
%%% Fundamentals Model
gubernatorial.nationalBiasIncumbency = -0.0226;
gubernatorial.nationalBiasMean = 0.0174;
gubernatorial.nationalBiasSigma = 0.0194;

gubernatorial.incumbency = 0.0603;
gubernatorial.presPrevSigma = 0.0713;
gubernatorial.presModelSigma = sqrt(0.0992^2 - 0.6043 * shared.districtSigma^2 - gubernatorial.presPrevSigma^2);
gubernatorial.prevModelSigma = sqrt(0.0881^2 - 0.6043 * shared.districtSigma^2 - gubernatorial.presPrevSigma^2);
gubernatorial.incumbencySigma = 0.0388;

%%% Polling
gubernatorial.nationalSigmaPoll = 0.0171;
gubernatorial.nationalQPoll = presidential.nationalQPoll * gubernatorial.nationalBiasSigma^2 / presidential.nationalSigma^2;

gubernatorial.stateSigmaPoll = sqrt(0.0276^2 - 0.6043 * shared.districtSigmaPoll^2);
gubernatorial.stateQPoll = gubernatorial.presModelSigma^2 / presidential.nationalSigma^2 * presidential.nationalQPoll;
gubernatorial.stateIncQPoll = gubernatorial.incumbencySigma^2 / presidential.nationalSigma^2 * presidential.nationalQPoll;

%% Save Data
save("Config.mat")