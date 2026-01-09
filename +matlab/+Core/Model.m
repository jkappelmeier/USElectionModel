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

            % Senate
            obj.loadSenatePrior();

            % Gubernatorial
            obj.loadGubernatorialPrior();

            % Cross-Race Correlations
            obj.loadSharedNationalCovariances();
            obj.loadSharedDistrictCovariances();
            obj.loadSharedIncumbencyCovariances();

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
                    H = double(electionType == electionInfoRedux.ElectionType' & election == electionInfoRedux.ElectionName' & geography == electionInfoRedux.GeographyName');

                    % Add bias factor
                    idxPres = election == "Presidential";
                    idxPresBias = electionInfoRedux.ElectionName == "Presidential" & electionInfoRedux.GeographyName == "National";
                    H(idxPres, idxPresBias) = 1;

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
            obj.covPoll = obj.covPoll/2 + obj.covPoll'/2;
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

            obj.covEst = obj.covEst/2 + obj.covEst'/2;
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
            districtCov = corr * obj.config.presidential.districtSigma^2;
            baseStates = regexprep(electionInfoPres.GeographyName, "-.*", "");
            stateCov = (baseStates == baseStates') * obj.config.presidential.stateSigma^2;
            stateCov(idxNat, idxNat) = 0;

            obj.covFund(idxPres, idxPres) = obj.covFund(idxPres,idxPres) + stateCov + districtCov;

            % Polling Uncertainty
            districtCovQPoll = corr * obj.config.presidential.districtQPoll;
            districtCovQBiasPoll = corr * obj.config.presidential.districtSigmaPoll^2;
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

            obj.xFund(idxGB) = xGBEst;
        end

        function obj = loadHousePrior(obj)
            arguments
                obj matlab.Core.Model
            end

            % Load Generic Ballot Fundamentals
            idxGB = obj.electionInfo.ElectionType == "Generic Ballot";
            inc = obj.electionInfo.IncumbencyFlag(idxGB);

            xHouseBiasEst = obj.config.house.nationalBiasMean + obj.config.house.nationalBiasIncumbency * inc;
            
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
            presHouseAdj = presDistrict - dot(presDistrict, totPresDistrict) / sum(totPresDistrict) + xHouseBiasEst + 0.5;
            xHouse(idxFlag) = presHouseAdj(idxFlag);

            houseVar = obj.config.house.presModelSigma^2;

            covHouse = zeros(N,N);
            covHouse(idxFlag, idxFlag) = eye(sum(idxFlag)) * houseVar;

            % Combine with incumbency
            xHouse(idxFlag) = xHouse(idxFlag) + electionInfoHouse.IncumbencyFlag(idxFlag) * obj.config.house.incumbency;

            obj.xFund(idxHouse) = xHouse;
            obj.covFund(idxHouse, idxHouse) = covHouse;


            % Polling
            obj.QBiasPoll(idxHouse, idxHouse) = eye(sum(idxHouse)) * obj.config.house.districtSigmaPoll^2;
            obj.QPoll(idxHouse, idxHouse) = eye(sum(idxHouse)) * obj.config.house.districtQPoll;

            % Add Previous House Election Data
            z = zeros(sum(idxHouse),1);
            measAvailFlag = true(sum(idxHouse),1);

            idxHousePrevBias = obj.preElectionData.ElectionType == "House" & obj.preElectionData.ElectionName == "Bias";
            idxHousePrevInc = obj.preElectionData.ElectionType == "House" & obj.preElectionData.ElectionName == "Incumbency";
            houseBias = obj.preElectionData.PreviousResult(idxHousePrevBias);
            houseInc = obj.preElectionData.PreviousResult(idxHousePrevInc);

            for i = 1:numel(z)
                incCur = electionInfoHouse.IncumbencyFlag(i);
                prevGeoName = obj.districtData.redistrictingData.OldDistrict(i);
                if incCur == 0 || idxFlag(i) == 0 || prevGeoName == "NaN"
                    measAvailFlag(i) = 0;
                else
                    idxPrev = obj.preElectionData.ElectionType == "House" & obj.preElectionData.GeographyName == prevGeoName;
                    zPrev = obj.preElectionData.PreviousResult(idxPrev);
                    if zPrev == 1 || zPrev == 0
                        measAvailFlag(i) = 0;
                    else
                        incPrev = obj.preElectionData.IncumbencyFlag(idxPrev);
                        adj = obj.districtData.redistrictingData.DistrictAdj(i);
                        z(i) = zPrev + (incCur * obj.config.house.incumbency - incPrev * houseInc) / 2 + adj + (xHouseBiasEst - houseBias) * 0.75;
                    end
                    
                end
            end

            z = z(measAvailFlag);

            % Create measurement noise matrix
            R = eye(numel(z)) * obj.config.house.prevModelSigma^2;

            % Do Measurment Update
            H = zeros(numel(z), numel(obj.xFund));
            idxMeas = idxHouse;
            idxMeas(idxHouse) = measAvailFlag;
            H(:, idxMeas) = eye(sum(idxMeas));

            y = z - H * obj.xFund;
            S = H * obj.covFund * H' + R;
            K = obj.covFund * H'/S;
            obj.xFund = obj.xFund + K * y;
            obj.covFund = (eye(numel(obj.xFund)) - K * H) * obj.covFund * (eye(numel(obj.xFund)) - K * H)' + K * R * K';

        end

        function obj = loadSenatePrior(obj)
            arguments
                obj matlab.Core.Model
            end

            % Load Generic Ballot Fundamentals
            idxGB = obj.electionInfo.ElectionType == "Generic Ballot";
            inc = obj.electionInfo.IncumbencyFlag(idxGB);

            xSenateBiasEst = obj.config.senate.nationalBiasMean + obj.config.senate.nationalBiasIncumbency * inc;

            % Load Senate Race Fundamentals
            idxSenate = obj.electionInfo.ElectionType == "Senate";
            electionInfoSenate = obj.electionInfo(idxSenate, :);

            noDemFlag = electionInfoSenate.CandidateD == "NaN";
            noRepFlag = electionInfoSenate.CandidateR == "NaN";
            idxFlag = ~noDemFlag & ~noRepFlag;
            obj.filterFlag(idxSenate) = idxFlag;

            N = sum(idxSenate);
            xSenate = zeros(N,1);

            idxPrevPres = obj.preElectionData.ElectionName == "Presidential";
            preElectionDataPres = obj.preElectionData(idxPrevPres, :);
            [~, idxState] = ismember(electionInfoSenate.GeographyName, preElectionDataPres.GeographyName);
            idxPrevPresNat = preElectionDataPres.GeographyName == "National";
            presState = preElectionDataPres.PreviousResult(idxState) - preElectionDataPres.PreviousResult(idxPrevPresNat) + 0.5 + xSenateBiasEst;
            xSenate(idxFlag) = presState(idxFlag);

            senateVar = obj.config.senate.presModelSigma^2;

            covSenate = zeros(N,N);
            covSenate(idxFlag, idxFlag) = eye(sum(idxFlag)) * senateVar;

            % Combine with incumbency
            xSenate(idxFlag) = xSenate(idxFlag) + electionInfoSenate.IncumbencyFlag(idxFlag) * obj.config.senate.incumbency;

            obj.xFund(idxSenate) = xSenate;
            obj.covFund(idxSenate, idxSenate) = covSenate;

            % Polling
            obj.QBiasPoll(idxSenate, idxSenate) = eye(sum(idxSenate)) * obj.config.senate.stateSigmaPoll^2;
            obj.QPoll(idxSenate, idxSenate) = eye(sum(idxSenate)) * obj.config.senate.stateQPoll;

            % Add Previous Senate Election Data
            z = zeros(sum(idxSenate),1);
            measAvailFlag = true(sum(idxSenate),1);

            idxSenatePrevBias = obj.preElectionData.ElectionType == "Senate" & obj.preElectionData.ElectionName == "Bias";
            idxSenatePrevInc = obj.preElectionData.ElectionType == "Senate" & obj.preElectionData.ElectionName == "Incumbency";
            senateBias = obj.preElectionData.PreviousResult(idxSenatePrevBias);
            senateInc = obj.preElectionData.PreviousResult(idxSenatePrevInc);
            idxSenatePrev = obj.preElectionData.ElectionType == "Senate" & obj.preElectionData.GeographyType == "State";
            preElectionDateSenate = obj.preElectionData(idxSenatePrev, :);

            for i = 1:numel(z)
                incCur = electionInfoSenate.IncumbencyFlag(i);
                if incCur == 0 || idxFlag(i) == 0
                    measAvailFlag(i) = 0;
                else
                    zPrev = preElectionDateSenate.PreviousResult(i);
                    if zPrev == 1 || zPrev == 0
                        measAvailFlag(i) = 0;
                    else
                        incPrev = preElectionDateSenate.IncumbencyFlag(i);
                        z(i) = zPrev + (incCur * obj.config.senate.incumbency - incPrev * senateInc) / 2 + (xSenateBiasEst - senateBias) * 0.75;
                    end
                    
                end
            end

            z = z(measAvailFlag);

            % Create measurement noise matrix
            R = eye(numel(z)) * obj.config.senate.prevModelSigma^2;

            % Do Measurment Update
            H = zeros(numel(z), numel(obj.xFund));
            idxMeas = idxSenate;
            idxMeas(idxSenate) = measAvailFlag;
            H(:, idxMeas) = eye(sum(idxMeas));

            y = z - H * obj.xFund;
            S = H * obj.covFund * H' + R;
            K = obj.covFund * H'/S;
            obj.xFund = obj.xFund + K * y;
            obj.covFund = (eye(numel(obj.xFund)) - K * H) * obj.covFund * (eye(numel(obj.xFund)) - K * H)' + K * R * K';

            obj.covFund(idxSenate, idxSenate) = obj.covFund(idxSenate, idxSenate) + eye(sum(idxSenate)) * obj.config.senate.presPrevSigma^2;
        end

        function obj = loadGubernatorialPrior(obj)
            arguments
                obj matlab.Core.Model
            end

            % Load Generic Ballot Fundamentals
            idxGB = obj.electionInfo.ElectionType == "Generic Ballot";
            inc = obj.electionInfo.IncumbencyFlag(idxGB);

            xGovBiasEst = obj.config.gubernatorial.nationalBiasMean + obj.config.gubernatorial.nationalBiasIncumbency * inc;

            % Load Gubernatorial Race Fundamentals
            idxGov = obj.electionInfo.ElectionType == "Gubernatorial";
            electionInfoGov = obj.electionInfo(idxGov, :);

            noDemFlag = electionInfoGov.CandidateD == "NaN";
            noRepFlag = electionInfoGov.CandidateR == "NaN";
            idxFlag = ~noDemFlag & ~noRepFlag;
            obj.filterFlag(idxGov) = idxFlag;

            N = sum(idxGov);
            xGov = zeros(N,1);

            idxPrevPres = obj.preElectionData.ElectionName == "Presidential";
            preElectionDataPres = obj.preElectionData(idxPrevPres, :);
            [~, idxState] = ismember(electionInfoGov.GeographyName, preElectionDataPres.GeographyName);
            idxPrevPresNat = preElectionDataPres.GeographyName == "National";
            presState = preElectionDataPres.PreviousResult(idxState) - preElectionDataPres.PreviousResult(idxPrevPresNat) + 0.5 + xGovBiasEst;
            xGov(idxFlag) = presState(idxFlag);

            govVar = obj.config.gubernatorial.presModelSigma^2;

            covGov = zeros(N,N);
            covGov(idxFlag, idxFlag) = eye(sum(idxFlag)) * govVar;

            % Combine with incumbency
            xGov(idxFlag) = xGov(idxFlag) + electionInfoGov.IncumbencyFlag(idxFlag) * obj.config.gubernatorial.incumbency;

            obj.xFund(idxGov) = xGov;
            obj.covFund(idxGov, idxGov) = covGov;

            % Polling
            obj.QBiasPoll(idxGov, idxGov) = eye(sum(idxGov)) * obj.config.gubernatorial.stateSigmaPoll^2;
            obj.QPoll(idxGov, idxGov) = eye(sum(idxGov)) * obj.config.gubernatorial.stateQPoll;

            % Add Previous Gubernatorial Election Data
            z = zeros(sum(idxGov),1);
            measAvailFlag = true(sum(idxGov),1);

            idxGovPrevBias = obj.preElectionData.ElectionType == "Gubernatorial" & obj.preElectionData.ElectionName == "Bias";
            idxGovPrevInc = obj.preElectionData.ElectionType == "Gubernatorial" & obj.preElectionData.ElectionName == "Incumbency";
            
            
            idxGovPrev = obj.preElectionData.ElectionType == "Gubernatorial" & obj.preElectionData.GeographyType == "State";
            preElectionDateGov = obj.preElectionData(idxGovPrev, :);

            prevPresYear = obj.preElectionData.ElectionYear(idxPrevPres);
            prevPresYear = prevPresYear(1);
            for i = 1:numel(z)
                incCur = electionInfoGov.IncumbencyFlag(i);
                if incCur == 0 || idxFlag(i) == 0
                    measAvailFlag(i) = 0;
                else
                    zPrev = preElectionDateGov.PreviousResult(i);
                    if zPrev == 1 || zPrev == 0
                        measAvailFlag(i) = 0;
                    else
                        prevYear = preElectionDateGov.ElectionYear(i);
                        incPrev = preElectionDateGov.IncumbencyFlag(i);
                        if prevYear > prevPresYear
                            govInc = obj.config.gubernatorial.incumbency;
                            govBias = xGovBiasEst;
                        else
                            govBias = obj.preElectionData.PreviousResult(idxGovPrevBias);
                            govInc = obj.preElectionData.PreviousResult(idxGovPrevInc);
                        end
                        
                        z(i) = zPrev + (incCur * obj.config.gubernatorial.incumbency - incPrev * govInc) / 2 + (xGovBiasEst - govBias);
                    end
                    
                end
            end

            z = z(measAvailFlag);

            % Create measurement noise matrix
            R = eye(numel(z)) * obj.config.gubernatorial.prevModelSigma^2;

            % Do Measurment Update
            H = zeros(numel(z), numel(obj.xFund));
            idxMeas = idxGov;
            idxMeas(idxGov) = measAvailFlag;
            H(:, idxMeas) = eye(sum(idxMeas));

            y = z - H * obj.xFund;
            S = H * obj.covFund * H' + R;
            K = obj.covFund * H'/S;
            obj.xFund = obj.xFund + K * y;
            obj.covFund = (eye(numel(obj.xFund)) - K * H) * obj.covFund * (eye(numel(obj.xFund)) - K * H)' + K * R * K';

            obj.covFund(idxGov, idxGov) = obj.covFund(idxGov, idxGov) + eye(sum(idxGov)) * obj.config.gubernatorial.presPrevSigma^2;
        end

        function obj = loadSharedNationalCovariances(obj)
            arguments
                obj matlab.Core.Model
            end

            % National Correlations
            idxPres = obj.electionInfo.ElectionType == "Presidential" & obj.electionInfo.GeographyType == "National";
            idxGB = obj.electionInfo.ElectionType == "Generic Ballot";
            idxHouse = obj.electionInfo.ElectionType == "House";
            idxSenate = obj.electionInfo.ElectionType == "Senate";
            idxGov = obj.electionInfo.ElectionType == "Gubernatorial";

            HBias = [idxPres, idxGB, idxHouse, idxSenate, idxGov];

            sigmaNatFund = [obj.config.presidential.nationalSigma; obj.config.genericBallot.nationalSigma; obj.config.house.nationalBiasSigma; obj.config.senate.nationalBiasSigma; obj.config.gubernatorial.nationalBiasSigma];
            rhoNationalFund = obj.config.shared.nationalFund;
            natBiasCov = sigmaNatFund * sigmaNatFund' .* rhoNationalFund;
            obj.covFund = obj.covFund + HBias * natBiasCov * HBias';

            sigmaNatPoll = [obj.config.presidential.nationalSigmaPoll; obj.config.genericBallot.nationalSigmaPoll; obj.config.house.nationalSigmaPoll; obj.config.senate.nationalSigmaPoll; obj.config.gubernatorial.nationalSigmaPoll];
            rhoNationalPoll = obj.config.shared.nationalPoll;
            natBiasPollCov = sigmaNatPoll * sigmaNatPoll' .* rhoNationalPoll;
            obj.QBiasPoll = obj.QBiasPoll + HBias * natBiasPollCov * HBias';

            sigmaNatQ = [obj.config.presidential.nationalQPoll; obj.config.genericBallot.nationalQPoll; obj.config.house.nationalQPoll; obj.config.senate.nationalQPoll; obj.config.gubernatorial.nationalQPoll].^0.5;
            natQPoll = sigmaNatQ * sigmaNatQ' .* rhoNationalFund;
            obj.QPoll = obj.QPoll + HBias * natQPoll * HBias';
        end

        function obj = loadSharedDistrictCovariances(obj)
            arguments
                obj matlab.Core.Model
            end

            % District Level Correlations
            rho = obj.districtData.corrMatrix;

            % Add in shared covariances
            covFundTemp = rho * obj.config.shared.districtSigma^2;
            QBiasPollTemp = rho * obj.config.shared.districtSigmaPoll^2;
            QPollTemp = rho * obj.config.shared.districtQPoll;

            % Combine shared covariances
            idxNonNat = obj.electionInfo.GeographyType ~= "National";
            geographyName = obj.electionInfo.GeographyName(idxNonNat);
            geographyType = obj.electionInfo.GeographyType(idxNonNat);
            H = obj.createMapping(geographyType, geographyName);

            obj.covFund(idxNonNat, idxNonNat) = obj.covFund(idxNonNat, idxNonNat) + H * covFundTemp * H';
            obj.QBiasPoll(idxNonNat, idxNonNat) = obj.QBiasPoll(idxNonNat, idxNonNat) + H * QBiasPollTemp * H';
            obj.QPoll(idxNonNat, idxNonNat) = obj.QPoll(idxNonNat, idxNonNat) + H * QPollTemp * H';
        end

        function loadSharedIncumbencyCovariances(obj)
            arguments
                obj matlab.Core.Model
            end

            % Combined shared covariances due to incumbency effects
            idxHouseInc = obj.electionInfo.ElectionType == "House" & obj.filterFlag;
            houseInc = 0 * obj.xFund;
            houseInc(idxHouseInc) = obj.electionInfo.IncumbencyFlag(idxHouseInc);
            houseIncSigma = obj.config.house.incumbencySigma;
            houseIncQ = obj.config.house.districtIncQPoll^0.5;

            idxSenateInc = obj.electionInfo.ElectionType == "Senate" & obj.filterFlag & obj.electionInfo.GeographyName ~= "National";
            senateInc = 0 * obj.xFund;
            senateInc(idxSenateInc) = obj.electionInfo.IncumbencyFlag(idxSenateInc);
            senateIncSigma = obj.config.senate.incumbencySigma;
            senateIncQ = obj.config.senate.stateIncQPoll^0.5;

            idxGovInc = obj.electionInfo.ElectionType == "Gubernatorial" & obj.filterFlag & obj.electionInfo.GeographyName ~= "National";
            govInc = 0 * obj.xFund;
            govInc(idxGovInc) = obj.electionInfo.IncumbencyFlag(idxGovInc);
            govIncSigma = obj.config.gubernatorial.incumbencySigma;
            govIncQ = obj.config.gubernatorial.stateIncQPoll^0.5;

            incRho = obj.config.shared.incCorr;
            HInc = [houseInc, senateInc, govInc];

            incSigma = [houseIncSigma; senateIncSigma; govIncSigma];
            incQ = [houseIncQ; senateIncQ; govIncQ];
            
            obj.covFund = obj.covFund + HInc * (incSigma * incSigma' .* incRho) * HInc';
            obj.QPoll = obj.QPoll + HInc * (incQ * incQ' .* incRho) * HInc';
            
        end
    end
end