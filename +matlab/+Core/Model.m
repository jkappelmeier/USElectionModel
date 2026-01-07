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

        end

        function obj = loadSharedNationalCovariances(obj)
            arguments
                obj matlab.Core.Model
            end

            % National Correlations
            idxPres = obj.electionInfo.ElectionType == "Presidential" & obj.electionInfo.GeographyType == "National";
            idxGB = obj.electionInfo.ElectionType == "Generic Ballot";
            idxSenate = obj.electionInfo.ElectionType == "Senate";
            idxGov = obj.electionInfo.ElectionType == "Gubernatorial";

            HBias = [idxPres, idxGB, idxSenate, idxGov];

            sigmaNatFund = [obj.config.presidential.nationalSigma; obj.config.genericBallot.nationalSigma; obj.config.senate.nationalBiasSigma; obj.config.gubernatorial.nationalBiasSigma];
            rhoNationalFund = obj.config.shared.nationalFund;
            natBiasCov = sigmaNatFund * sigmaNatFund' .* rhoNationalFund;
            obj.covFund = obj.covFund + HBias * natBiasCov * HBias';

            sigmaNatPoll = [obj.config.presidential.nationalSigmaPoll; obj.config.genericBallot.nationalSigmaPoll; obj.config.senate.nationalSigmaPoll; obj.config.gubernatorial.nationalSigmaPoll];
            rhoNationalPoll = obj.config.shared.nationalPoll;
            natBiasPollCov = sigmaNatPoll * sigmaNatPoll' .* rhoNationalPoll;
            obj.QBiasPoll = obj.QBiasPoll + HBias * natBiasPollCov * HBias';

            sigmaNatQ = [obj.config.presidential.nationalQPoll; obj.config.genericBallot.nationalQPoll; obj.config.senate.nationalQPoll; obj.config.gubernatorial.nationalQPoll].^0.5;
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
            houseInc = obj.electionInfo.IncumbencyFlag(idxHouseInc);
            houseIncSigma = obj.config.house.incumbencySigma;
            houseIncSigmaVec = houseIncSigma * houseInc;
            houseIncSigmaVecQ = obj.config.house.districtIncQPoll^0.5 * houseInc;

            idxSenateInc = obj.electionInfo.ElectionType == "Senate" & obj.filterFlag & obj.electionInfo.GeographyName ~= "National";
            senateInc = obj.electionInfo.IncumbencyFlag(idxSenateInc);
            senateIncSigma = obj.config.senate.incumbencySigma;
            senateIncSigmaVec = senateIncSigma * senateInc;
            senateIncSigmaVecQ = obj.config.senate.stateIncQPoll^0.5 * senateInc;

            idxGovInc = obj.electionInfo.ElectionType == "Gubernatorial" & obj.filterFlag & obj.electionInfo.GeographyName ~= "National";
            govInc = obj.electionInfo.IncumbencyFlag(idxGovInc);
            govIncSigma = obj.config.gubernatorial.incumbencySigma;
            govIncSigmaVec = govIncSigma * govInc;
            govIncSigmaVecQ = obj.config.gubernatorial.stateIncQPoll^0.5 * govInc;

            incRho = obj.config.shared.incCorr;

            idx  = {idxHouseInc, idxSenateInc, idxGovInc};
            sig  = {houseIncSigmaVec(:), senateIncSigmaVec(:), govIncSigmaVec};
            
            for i = 1:numel(idx)
                for j = i:numel(idx)
                    blk = (sig{i} * sig{j}.') * incRho(i,j);
                    obj.covFund(idx{i}, idx{j}) = obj.covFund(idx{i}, idx{j}) + blk;
                    if j>i, obj.covFund(idx{j}, idx{i}) = obj.covFund(idx{j}, idx{i}) + blk.'; end
                end
            end

            sig  = {houseIncSigmaVecQ(:), senateIncSigmaVecQ(:), govIncSigmaVecQ};
            
            for i = 1:numel(idx)
                for j = i:numel(idx)
                    blk = (sig{i} * sig{j}.') * incRho(i,j);
                    obj.QPoll(idx{i}, idx{j}) = obj.QPoll(idx{i}, idx{j}) + blk;
                    if j>i, obj.QPoll(idx{j}, idx{i}) = obj.QPoll(idx{j}, idx{i}) + blk.'; end
                end
            end
            
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