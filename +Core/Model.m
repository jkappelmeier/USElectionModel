classdef Model < handle
    properties (Access = public)
        name (1,1) string
        electionInfo table
        preElectionData table
        correlationData struct
        config struct

        xFund (:,1) double
        covFund (:,:) double

        hDistrict2State (:,:) double
        filterFlag (:,1) logical
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

            % ELECTION AND PRE-ELECTION DATE MUST HAVE DSAME ELECTION AND
            % GEOGRAPHY

            N = size(obj.electionInfo, 1);
            obj.xFund = zeros(N,1);
            obj.covFund = zeros(N,N);
            obj.filterFlag = ones(N,1);

            % Presidential
            % Estimates
            idxNat = obj.electionInfo.GeographyType == 'National';
            incPres = obj.electionInfo.IncumbencyFlag(idxNat);
            natEst = obj.config.nationalMean + obj.config.incumbency * incPres;

            prevNat = obj.preElectionData.PreviousResult(idxNat);

            obj.xFund(~idxNat) = obj.preElectionData.PreviousResult(~idxNat) - prevNat;
            obj.xFund(idxNat) = natEst;

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
                obj Core.Model
                nSims (1,1) int32
            end

            xSims = mvnrnd(obj.xFund, obj.covFund, nSims)';
        end
    end
end