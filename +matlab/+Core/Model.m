classdef Model < handle
    properties (Access = public)
        name (1,1) string
        electionInfo table
        preElectionData table
        districtData struct
        config struct

        time (1,:) double
        xFund (:,1) double
        covFund (:,:) double

        filterFlag (:,1) logical

        QPoll (:,:) double
        QBiasPoll (:,:) double

        pollTable table;
        xPoll (:,1) double
        covPoll (:,:) double
        availFlagPoll (:,:) logical

        xEst (:,1) double
        covEst (:,:) double
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
                districtData, config)
            arguments
                name (1,1) string
                electionInfo (1,1) string
                preElectionData (1,1) string
                districtData (1,1) string
                config (1,1) string
            end

            % Load Data
            obj.name = name;
            obj.electionInfo = readtable(electionInfo, "TextType", "string");
            obj.districtData = load(districtData);
            obj.config = load(config);
            obj.preElectionData = readtable(preElectionData, "TextType", "string");

            % Pre-Allocate States
            N = size(obj.electionInfo, 1);
            obj.xFund = zeros(N,1);
            obj.covFund = sparse(zeros(N,N));
            obj.filterFlag = ones(N,1);
            
            t0 = days(obj.config.electionDate - obj.config.startDate);
            obj.time = t0:-1:0;
            obj.availFlagPoll = zeros(N,length(obj.time));

            % Polling Error Initialization
            obj.QPoll = sparse(zeros(N,N));
            obj.QBiasPoll = sparse(zeros(N,N));

            % Presidential
            obj.loadPresidentialPrior();

            % Generic Ballot
            obj.loadGenericBallotPrior();

            % House
            obj.loadHousePrior();

            % Cross-Race Correlations
            obj.loadSharedCovariances();

            % Eliminate Uncessary Rows
            obj.QPoll = obj.QPoll(obj.filterFlag, obj.filterFlag);
            obj.QBiasPoll = obj.QBiasPoll(obj.filterFlag, obj.filterFlag);

            % Ensure symmetric
            obj.covFund = (obj.covFund + obj.covFund') / 2;
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

            xSimsTemp = mvnrnd(obj.xEst(obj.filterFlag), obj.covEst(obj.filterFlag, obj.filterFlag), nSims)';
            xSims = repmat(obj.xEst, 1, nSims);
            xSims(obj.filterFlag, :) = xSimsTemp;

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

            obj.pollTable = table(pollData.ElectionType, pollData.Election, pollData.Geography, pollData.GeographyType, ...
                t, z, R, 'VariableNames', {'ElectionType', 'Election', 'Geography', 'GeographyType', ...
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
            x0 = obj.xFund(obj.filterFlag);
            cov0 = obj.covFund(obj.filterFlag, obj.filterFlag) * 1000000;

            xPollCur = x0;
            covPollCur = cov0;

            N = length(obj.time);
            n = length(x0);
            obj.xPoll = zeros(n,1);
            obj.covPoll = zeros(n,n);
            
            % Some Preallocation for time saving
            electionInfoRedux = obj.electionInfo(obj.filterFlag, :);
            availFlagRedux = obj.availFlagPoll(obj.filterFlag, :);

            % Loop through times
            for i = 1:N
                % Get poll results at current time
                idx = obj.pollTable.tTillElection == obj.time(i);
                
                % If current time has polls
                if sum(idx) > 0
                    zVec = obj.pollTable.zVec(idx);
                    R = diag(obj.pollTable.R(idx));
                    electionType = obj.pollTable.ElectionType(idx);
                    election = obj.pollTable.Election(idx);
                    geography = obj.pollTable.Geography(idx);

                    % Generate H matrix for current time step
                    H = double(election == electionInfoRedux.ElectionName' & geography == electionInfoRedux.GeographyName');

                    % Add bias factor
                    idxPres = election == "Presidential";
                    idxPresBias = electionInfoRedux.ElectionName == "Presidential" & electionInfoRedux.GeographyName == "National";
                    H(idxPres, idxPresBias) = 1;

                    idxHouse = electionType == "House";
                    idxGB = electionInfoRedux.ElectionName == "Generic Ballot";
                    H(idxHouse, idxGB) = 1;

                    availFlagRedux(:,i) = any(H > 0, 1)';

                    % Update
                    yVec = zVec - H * xPollCur;
                    K = covPollCur * H' / (H*covPollCur*H' + R);
                    xPollCur = xPollCur + K*yVec;
                    covPollCur = (eye(n) - K*H) * covPollCur * (eye(n) - K*H)' + K*R*K'; 
                end

                covPollCur = covPollCur + obj.QPoll;
            end

            obj.covPoll = covPollCur + obj.QBiasPoll;
            obj.xPoll = xPollCur;

            obj.availFlagPoll(obj.filterFlag, :) = cummax(availFlagRedux, 2);
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

            N = length(x);

            pollAvailFlag = obj.availFlagPoll(obj.filterFlag,end);
            stateAvailFlag = obj.availFlagPoll(:,end) & obj.filterFlag;
            if any(stateAvailFlag)
                zCur = z(pollAvailFlag);
                RCur = R(pollAvailFlag,pollAvailFlag);
                H = eye(N);
                H(~stateAvailFlag, :) = [];

                y = zCur - H*x;
                K = P*H'/(H*P*H'+RCur);

                obj.xEst = x + K * y;
                obj.covEst = (eye(N)-K*H)*P*(eye(N)-K*H)'+K*RCur*K';
            else
                obj.xEst = x;
                obj.covEst = P;
            end
        end
    end

    methods (Access = private)

        % Create the prior for the Presidential Model
        %
        % Input:
        %   obj - Instance of this object
        % Output:
        %   obj - Instance of this object
        function obj = loadPresidentialPrior(obj)
            arguments
                obj matlab.Core.Model
            end
            
            idxPres = obj.electionInfo.ElectionType == "Presidential";
            electionInfoPres = obj.electionInfo(idxPres,:);
            idxNat = obj.electionInfo.GeographyType == "National" & obj.electionInfo.ElectionType == "Presidential";
            incPres = obj.electionInfo.IncumbencyFlag(idxNat);
            natEst = obj.config.presidential.nationalMean + obj.config.presidential.incumbency * incPres;

            idxPrevNat = obj.preElectionData.GeographyType == "National" & obj.preElectionData.ElectionType == "Presidential";
            prevNat = obj.preElectionData.PreviousResult(idxPrevNat);
            idxPrevPres = obj.preElectionData.ElectionType == "Presidential";

            obj.xFund(idxNat) = natEst;
            obj.xFund(~idxNat & idxPres) = obj.preElectionData.PreviousResult(~idxPrevNat & idxPrevPres) - prevNat;

            % Uncertainty
            hDistrict2State = obj.createMapping(electionInfoPres.GeographyType, electionInfoPres.GeographyName);

            % Initial Uncertainty
            corr = hDistrict2State * obj.districtData.corrMatrix * hDistrict2State';
            obj.covFund(idxNat, idxNat) = obj.config.presidential.nationalSigma^2;
            districtCov = corr * (obj.config.presidential.districtSigma^2 - obj.config.shared.districtSigma^2);
            baseStates = regexprep(electionInfoPres.GeographyName, "-.*", "");
            stateCov = (baseStates == baseStates') * obj.config.presidential.stateSigma^2;
            stateCov(idxNat, idxNat) = 0;

            obj.covFund(idxPres, idxPres) = obj.covFund(idxPres,idxPres) + stateCov + districtCov;

            % Polling Uncertainty
            obj.QPoll(idxNat, idxNat) = obj.config.presidential.nationalQPoll;
            obj.QBiasPoll(idxNat, idxNat) = obj.config.presidential.nationalSigmaPoll^2;
            districtCovQPoll = corr * (obj.config.presidential.districtQPoll - obj.config.shared.districtQPoll);
            districtCovQBiasPoll = corr * (obj.config.presidential.districtSigmaPoll^2 - obj.config.shared.districtSigmaPoll^2);
            stateCovQPoll = (baseStates == baseStates') * obj.config.presidential.stateQPoll;
            stateCovQBiasPoll = (baseStates == baseStates') * obj.config.presidential.stateSigmaPoll^2;
            stateCovQPoll(idxNat, idxNat) = 0;
            stateCovQBiasPoll(idxNat, idxNat) = 0;

            obj.QPoll(idxPres, idxPres) = obj.QPoll(idxPres, idxPres) + stateCovQPoll + districtCovQPoll;
            obj.QBiasPoll(idxPres, idxPres) = obj.QBiasPoll(idxPres, idxPres) + stateCovQBiasPoll + districtCovQBiasPoll;

        end

        function obj = loadGenericBallotPrior(obj)
            arguments
                obj matlab.Core.Model
            end
            
            % Load Generic Ballot Fundamentals
            idxGB = obj.electionInfo.ElectionType == "Generic Ballot";
            inc = obj.electionInfo.IncumbencyFlag(idxGB);

            xGBEst = obj.config.genericBallot.nationalMean + obj.config.genericBallot.nationalIncumbency * inc;
            covGBEst = obj.config.genericBallot.nationalSigma^2;

            obj.xFund(idxGB) = xGBEst;
            obj.covFund(idxGB,idxGB) = covGBEst;

            % Load Generic Ballot Polling
            obj.QBiasPoll(idxGB,idxGB) = obj.config.genericBallot.nationalSigmaPoll^2;
            obj.QPoll(idxGB, idxGB) = obj.config.genericBallot.nationalQPoll;

        end

        function obj = loadHousePrior(obj)
            arguments
                obj matlab.Core.Model
            end
            
            % Load House Fundamentals
            idxHouse = obj.electionInfo.ElectionType == "House";
            electionInfoHouse = obj.electionInfo(idxHouse, :);

            noDemFlag = electionInfoHouse.CandidateD == "NaN";
            noRepFlag = electionInfoHouse.CandidateR == "NaN";
            idxFlag = ~noDemFlag & ~noRepFlag;
            obj.filterFlag(idxHouse) = idxFlag;

            N = sum(idxHouse);
            xHouse = zeros(N,1);

            H = obj.createMapping(electionInfoHouse.GeographyType, electionInfoHouse.GeographyName);
            presDistrict = H * obj.districtData.districtInfo.PresResult;
            totPresDistrict = H * obj.districtData.districtInfo.TotalVote;
            presHouseAdj = presDistrict - dot(presDistrict, totPresDistrict) / sum(totPresDistrict);
            xHouse(idxFlag) = presHouseAdj(idxFlag);

            houseVar = obj.config.house.presModelSigma^2 - obj.config.shared.districtSigma^2;

            covHouse = zeros(N,N);
            covHouse(idxFlag, idxFlag) = eye(sum(idxFlag)) * houseVar;

            % Combine with incumbency
            xHouse(idxFlag) = xHouse(idxFlag) + electionInfoHouse.IncumbencyFlag(idxFlag) * obj.config.house.incumbency;

            obj.xFund(idxHouse) = xHouse;
            obj.covFund(idxHouse, idxHouse) = covHouse;


            % Polling
            obj.QBiasPoll(idxHouse, idxHouse) = eye(sum(idxHouse)) * (obj.config.house.districtSigmaPoll^2 - obj.config.shared.districtSigmaPoll^2);
            obj.QPoll(idxHouse, idxHouse) = eye(sum(idxHouse)) * (obj.config.house.districtQPoll - obj.config.shared.districtQPoll);

        end

        function obj = loadSharedCovariances(obj)
            arguments
                obj matlab.Core.Model
            end

            % National Correlations
            idxPres = obj.electionInfo.ElectionType == "Presidential" & obj.electionInfo.GeographyType == "National";
            idxGB = obj.electionInfo.ElectionType == "Generic Ballot";
            idx = idxPres | idxGB;

            sigmaNatPoll = [obj.config.presidential.nationalSigmaPoll; obj.config.genericBallot.nationalSigmaPoll];
            rhoNationalPoll = obj.config.shared.nationalPoll;
            obj.QBiasPoll(idx,idx) = sigmaNatPoll * sigmaNatPoll' .* rhoNationalPoll;

            % District Level Correlations
            rho = obj.districtData.corrMatrix;

            % Add in shared covariances
            covFundTemp = rho * obj.config.shared.districtSigma^2;
            QBiasPollTemp = rho * obj.config.shared.districtSigmaPoll^2;
            QPollTemp = rho * obj.config.shared.districtQPoll;

            % Combine shared covariances
            idxNonIncNat = obj.electionInfo.GeographyType ~= "National" & obj.electionInfo.ElectionName ~= "Incumbency" & obj.filterFlag;
            geographyName = obj.electionInfo.GeographyName(idxNonIncNat);
            geographyType = obj.electionInfo.GeographyType(idxNonIncNat);
            H = obj.createMapping(geographyType, geographyName);

            obj.covFund(idxNonIncNat, idxNonIncNat) = obj.covFund(idxNonIncNat, idxNonIncNat) + H * covFundTemp * H';
            obj.QBiasPoll(idxNonIncNat, idxNonIncNat) = obj.QBiasPoll(idxNonIncNat, idxNonIncNat) + H * QBiasPollTemp * H';
            obj.QPoll(idxNonIncNat, idxNonIncNat) = obj.QPoll(idxNonIncNat, idxNonIncNat) + H * QPollTemp * H';

            % Combined shared covariances due to incumbency effects
            idxInc = obj.electionInfo.ElectionType == "House" & obj.filterFlag;
            houseInc = obj.electionInfo.IncumbencyFlag(idxInc);
            houseIncSigma = obj.config.house.incumbencySigma;
            houseIncSigmaVec = houseIncSigma * houseInc;
            houseIncSigmaVecQ = obj.config.house.districtIncQPoll^0.5 * houseInc;

            obj.covFund(idxInc, idxInc) = obj.covFund(idxInc, idxInc) + houseIncSigmaVec * houseIncSigmaVec';
            obj.QPoll(idxInc, idxInc) = obj.QPoll(idxInc, idxInc) + houseIncSigmaVecQ * houseIncSigmaVecQ';
            
            
        end

        function H = createMapping(obj, geographyTypes, geographyNames)
            arguments
                obj matlab.Core.Model
                geographyTypes (:,1) string
                geographyNames (:,1) string
            end

            m = length(geographyTypes);
            N = length(obj.districtData.districtInfo.District);
            districtNames = obj.districtData.districtInfo.District;
            stateNames = extractBefore(districtNames, "-");
            totalVote = obj.districtData.districtInfo.TotalVote;

            H = zeros(m,N);

            for i = 1:m
                if geographyTypes(i) == "Congressional District"
                    idx = districtNames == geographyNames(i);
                    H(i, idx) = 1;
                elseif geographyTypes(i) == "State"
                    idx = strcmp(stateNames, geographyNames(i));
                    prop = totalVote(idx)/sum(totalVote(idx));
                    H(i, idx) = prop';
                end
            end
        end
    end
end