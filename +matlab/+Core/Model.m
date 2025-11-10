classdef Model < handle
    properties (Access = public)
        name (1,1) string
        electionInfo table
        preElectionData table
        correlationData struct
        config struct

        time (1,:) double
        xFund (:,1) double
        covFund (:,:) double

        hDistrict2State (:,:) double
        filterFlag (:,1) logical

        QPoll (:,:) double
        QBiasPoll (:,:) double

        pollTable table;
        xPoll (:,:) double
        covPoll (:,:,:) double
        availFlagPoll (:,:) logical

        xEst (:,:) double
        covEst (:,:,:) double
    end

    methods (Access = public)
        % Constructor which creates xFund and covFund
        %
        % Inputs:
        %   name            - name of model
        %   electionInfo    - path to election info csv file
        %   preElectionData - path to pre election data csv file
        %   correlationData - path to correlation data mat file
        %   config          - path to config for the election, mat file
        %
        % Outputs:
        function obj = Model(name, electionInfo, preElectionData, ...
                correlationData, config)
            arguments
                name (1,1) string
                electionInfo (1,1) string
                preElectionData (1,1) string
                correlationData (1,1) string
                config (1,1) string
            end

            % Load Data
            obj.name = name;
            obj.electionInfo = readtable(electionInfo, "TextType", "string");
            obj.correlationData = load(correlationData);
            obj.config = load(config);
            obj.preElectionData = readtable(preElectionData, "TextType", "string");

            % ELECTION AND PRE-ELECTION DATE MUST HAVE SAME ELECTION AND
            % GEOGRAPHY

            N = size(obj.electionInfo, 1);
            obj.xFund = zeros(N,1);
            obj.covFund = zeros(N,N);
            obj.filterFlag = ones(N,1);

            % Presidential
            t0 = days(obj.config.electionDate - obj.config.startDate);
            obj.time = t0:-1:0;
            % Estimates
            idxNat = obj.electionInfo.GeographyType == 'National';
            incPres = obj.electionInfo.IncumbencyFlag(idxNat);
            natEst = obj.config.nationalMean + obj.config.incumbency * incPres;

            prevNat = obj.preElectionData.PreviousResult(idxNat);

            obj.xFund(idxNat) = natEst;
            obj.xFund(~idxNat) = obj.preElectionData.PreviousResult(~idxNat) - prevNat;

            % Uncertainty

            districtNames = obj.correlationData.corrInfo.District;
            stateNames = extractBefore(districtNames, "-");
            totalVote = obj.correlationData.corrInfo.TotalVote;
            nDistricts = numel(districtNames);
            obj.hDistrict2State = zeros(N, nDistricts);
            
            for i = 1:N
                if obj.electionInfo.GeographyType(i) == "Congressional District"
                    idx = districtNames == obj.electionInfo.GeographyName(i);
                    obj.hDistrict2State(i, idx) = 1;
                elseif obj.electionInfo.GeographyType(i) == "State"
                    idx = strcmp(stateNames, obj.electionInfo.GeographyName(i));
                    prop = totalVote(idx)/sum(totalVote(idx));
                    obj.hDistrict2State(i, idx) = prop';
                end
            end

            corr = obj.hDistrict2State * obj.correlationData.corrMatrix * obj.hDistrict2State';
            obj.covFund(idxNat, idxNat) = obj.config.nationalSigma^2;
            districtCov = corr * obj.config.districtSigma^2;
            baseStates = regexprep(obj.electionInfo.GeographyName, "-.*", "");
            stateCov = (baseStates == baseStates') * obj.config.stateSigma^2;
            stateCov(idxNat, idxNat) = 0;

            obj.covFund = obj.covFund + stateCov + districtCov;
            
            % Ensure symmetric
            obj.covFund = (obj.covFund + obj.covFund') / 2;

            % Polling Error Initialization
            obj.QPoll = zeros(N,N);
            obj.QBiasPoll = zeros(N,N);

            obj.QPoll(idxNat, idxNat) = obj.config.nationalQPoll;
            obj.QBiasPoll(idxNat, idxNat) = obj.config.nationalSigmaPoll^2;
            districtCovQPoll = corr * obj.config.districtQPoll;
            districtCovQBiasPoll = corr * obj.config.districtSigmaPoll^2;
            stateCovQPoll = (baseStates == baseStates') * obj.config.stateQPoll;
            stateCovQBiasPoll = (baseStates == baseStates') * obj.config.stateSigmaPoll^2;
            stateCovQPoll(idxNat, idxNat) = 0;
            stateCovQBiasPoll(idxNat, idxNat) = 0;

            obj.QPoll = obj.QPoll + stateCovQPoll + districtCovQPoll;
            obj.QBiasPoll = obj.QBiasPoll + stateCovQBiasPoll + districtCovQBiasPoll;

            obj.QPoll = (obj.QPoll + obj.QPoll') / 2;
            obj.QBiasPoll = (obj.QBiasPoll + obj.QBiasPoll') / 2;
        end

        % Simulate election using the state vector and covariance
        %
        % Inputs:
        %   obj     - Model objet
        %   nSims   - number of sims to generate
        %
        % Output:
        %   xSims   - Matrix with simulation results
        function xSims = simulate(obj, nSims)
            arguments
                obj matlab.Core.Model
                nSims (1,1) int32
            end

            xSims = mvnrnd(obj.xEst(:,end), obj.covEst(:,:,end), nSims)';
        end

        % Add polls to model
        %
        % Inputs:
        %   obj         - Model object
        %   pollFile    - csv file containing polls
        %
        % Output:
        %   obj - Model object
        function addPolls(obj, pollFile)
            arguments
                obj matlab.Core.Model
                pollFile (1,1) string
            end

            tElection = obj.config.electionDate;

            pollData = readtable(pollFile, "TextType", "string");
            N = height(pollData);

            t = zeros(N,1);
            z = zeros(N,1);
            R = zeros(N,1);
            for i = 1:N
                tStart = days(tElection - pollData.DateStart(i));
                tEnd = days(tElection - pollData.DateEnd(i));
                tMid = floor(mean([tStart, tEnd]));
                t(i) = tMid;

                dem = str2double(erase(pollData.Democrat(i), "%"))/100;
                rep = str2double(erase(pollData.Republican(i), "%"))/100;

                z(i) = dem / (dem + rep);
                R(i) = 0.25/(pollData.SampleSize(i)*(dem+rep));
            end

            obj.pollTable = table(pollData.Election, pollData.Geography, ...
                t, z, R, 'VariableNames', {'Election', 'Geography', ...
                'tTillElection', 'zVec', 'R'});

            obj.pollTable = sortrows(obj.pollTable, ...
                "tTillElection", "descend");
        end

        % Run Polling Average
        %
        % Input
        %   obj - Model object
        %
        % Output
        %   obj - Model Object
        function obj = runPollingAverage(obj)
            arguments
                obj matlab.Core.Model
            end

            % Initial state and cov
            x0 = obj.xFund;
            cov0 = obj.covFund * 1000000;

            xPollCur = x0;
            covPollCur = cov0;

            N = length(obj.time);
            n = length(x0);
            obj.xPoll = zeros(n,N);
            obj.covPoll = zeros(n,n,N);
            obj.availFlagPoll = zeros(n,N);

            % Loop through times
            for i = 1:N
                % Get poll results at current time
                idx = obj.pollTable.tTillElection == obj.time(i);
                
                % If current time has polls
                if sum(idx) > 0
                    zVec = obj.pollTable.zVec(idx);
                    R = diag(obj.pollTable.R(idx));
                    election = obj.pollTable.Election(idx);
                    geography = obj.pollTable.Geography(idx);

                    % Generate H matrix for current time step
                    H = double(election == obj.electionInfo.ElectionName' & geography == obj.electionInfo.GeographyName');

                    % Add bias factor
                    idxPres = election == "Presidential";
                    idxPresBias = obj.electionInfo.ElectionName == "Presidential" & obj.electionInfo.GeographyName == "National";
                    H(idxPres, idxPresBias) = 1;

                    obj.availFlagPoll(:,i) = any(H, 1)';

                    % Update
                    yVec = zVec - H * xPollCur;
                    K = covPollCur * H' / (H*covPollCur*H' + R);
                    xPollCur = xPollCur + K*yVec;
                    covPollCur = (eye(n) - K*H) * covPollCur * (eye(n) - K*H)' + K*R*K';
                else
                    
                end

                obj.xPoll(:,i) = xPollCur;
                obj.covPoll(:,:,i) = covPollCur + obj.QPoll * obj.time(i) + obj.QBiasPoll;

                covPollCur = covPollCur + obj.QPoll;
            end

            obj.availFlagPoll = cummax(obj.availFlagPoll, 2);
        end

        % Combine the polls with the fundamentals
        %
        % Input:
        %   obj - Model object
        %
        % Output:
        %   obj - Model object
        function obj = generateEstimate(obj)
            arguments
                obj matlab.Core.Model
            end

            x = obj.xFund;
            P = obj.covFund;

            z = obj.xPoll;
            R = obj.covPoll;

            nTimes = length(obj.time);
            N = length(x);
            
            for i = 1:nTimes
                availFlagCur = obj.availFlagPoll(:,i);
                if any(availFlagCur)
                    zCur = z(availFlagCur,i);
                    RCur = R(availFlagCur,availFlagCur,i);
                    H = eye(N);
                    H(~availFlagCur, :) = [];
    
                    y = zCur - H*x;
                    K = P*H'/(H*P*H'+RCur);
    
                    obj.xEst(:,i) = x + K * y;
                    obj.covEst(:,:,i) = (eye(N)-K*H)*P*(eye(N)-K*H)'+K*RCur*K';
                else
                    obj.xEst(:,i) = x;
                    obj.covEst(:,:,i) = P;
                end
            end
        end
    end
end