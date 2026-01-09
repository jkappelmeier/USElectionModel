classdef VisualizeResults
    methods (Static)
        function printResults(model, xSims)
            arguments
                model matlab.Core.Model
                xSims (:,:) double
            end

            ecVotes = matlab.Visualization.VisualizeResults.printPresidentialResults(model, xSims);
            matlab.Visualization.VisualizeResults.printHouseResults(model, xSims);
            matlab.Visualization.VisualizeResults.printSenateResults(model, xSims, ecVotes);
            matlab.Visualization.VisualizeResults.printGubernatorialResults(model);
        end

        function ecVotes = printPresidentialResults(model, xSims)
            arguments
                model matlab.Core.Model
                xSims (:,:) double
            end

            % Use only presidential Results
            presIdx = model.electionInfo.ElectionType == "Presidential";
            electionInfo = model.electionInfo(presIdx,:);
            xEstPres = model.xEst(presIdx,:);
            covEstPres = model.covEst(presIdx,presIdx,:);
            preElectionDataPres = model.preElectionData(presIdx,:);
            xSimsPres = xSims(presIdx,:);


            % Get transformation from x to y
            idxNat = electionInfo.GeographyType == "National";
            N = numel(xEstPres);
            hX2Y = zeros(N, N - 1);
            hX2Y(~idxNat, :) = eye(N-1);
            hX2Y(idxNat, :) = ones(1, N-1);

            xEst = hX2Y' * xEstPres;
            xSimsPres = hX2Y' * xSimsPres;
            covEst = hX2Y' * covEstPres * hX2Y;
            score = electionInfo.Score(~idxNat);
            geographyType = electionInfo.GeographyType(~idxNat);

            % Get electoral votes per run
            nSims = size(xSimsPres, 2);
            ecVotes = score' * (xSimsPres > 0.5);
            ecMean = mean(ecVotes);
            presECWin = sum(ecVotes >= 270) / nSims;
            presECLose = sum(ecVotes < 268) / nSims;
            presECTie = 1 - presECWin - presECLose;

            % Get popular vote
            hState2Popular = 0 * xEst;
            idxState = geographyType == "State";
            prevVoteTot = preElectionDataPres.PreviousVote(~idxNat);
            hState2Popular(idxState) = prevVoteTot(idxState);
            hState2Popular = hState2Popular / sum(hState2Popular);
            meanPopVote = xEst' * hState2Popular;
            
            sigmaPopVote = sqrt(hState2Popular' * covEst * ...
                hState2Popular);

            % Display the results
            dCand = electionInfo.CandidateD(idxNat);
            rCand = electionInfo.CandidateR(idxNat);

            disp("--------------------------------------------------------")
            disp("Presidential Election:")
            disp("  " + dCand + ": " + num2str(ecMean,"%.2f") + " Electoral Votes | Chance of Winning: " + num2str(presECWin*100, "%.2f") + "%");
            disp("  " + rCand + ": " + num2str(538-ecMean,"%.2f") + " Electoral Votes | Chance of Winning: " + num2str(presECLose*100, "%.2f") + "%");
            disp("  Chance of Tie: " + num2str(presECTie*100,"%.2f") + "%");
            disp(" ")

            disp("Popular Vote:")
            popVoteChance = normcdf(meanPopVote, 0.5, sigmaPopVote);
            disp("  " + dCand + ": " + num2str(meanPopVote*100,"%.2f") + "% | Chance of Winning: " + num2str(popVoteChance * 100,"%.2f") + "%");
            disp("  " + rCand + ": " + num2str(100-meanPopVote*100,"%.2f") + "% | Chance of Winning: " + num2str(100-popVoteChance * 100,"%.2f") + "%");
            disp("--------------------------------------------------------")

            geographyNames = electionInfo.GeographyName(~idxNat);
            for i = 1:N-1
                xEstCur = xEst(i);
                sigmaEst = sqrt(covEst(i,i));
                chance = normcdf(xEstCur, 0.5, sigmaEst);
                geoName = geographyNames(i);
                if geographyType(i) == "Congressional District"
                    state = extractBefore(geoName, "-");
                    num = extractAfter(geoName, "-");
                    firstDig = str2double(extractBetween(num, 1, 1));
                    lastDig = str2double(extractBetween(num, 2, 2));
                    if lastDig == 1 && firstDig ~=1
                        numStr = string(str2double(num)) + "st";
                    elseif lastDig == 2 && firstDig ~=1
                        numStr = string(str2double(num)) + "nd";
                    elseif lastDig == 3 && firstDig ~=1
                        numStr = string(str2double(num)) + "rd";
                    else
                        numStr = string(str2double(num)) + "th";
                    end

                    geoName = state + " " + numStr;
                end
                disp(geoName + " (" + num2str(score(i),"%.0f") + " Electoral Votes):")
                disp("  " + dCand + ": " + num2str(xEstCur*100,"%.2f") + "% | Chance of Winning: " + num2str(chance*100,"%.2f") + "%")
                disp("  " + rCand + ": " + num2str((1-xEstCur)*100,"%.2f") + "% | Chance of Winning: " + num2str((1-chance)*100,"%.2f") + "%")
                disp(" ")
            end
        end

        function printHouseResults(model, xSims)
            arguments
                model matlab.Core.Model
                xSims (:,:) double
            end

            % Model House Results
            houseIdx = model.electionInfo.ElectionType == "House";
            electionInfo = model.electionInfo(houseIdx,:);
            xEstHouse = model.xEst(houseIdx,:);
            covEstHouse = model.covEst(houseIdx,houseIdx);
            xSimsHouse = xSims(houseIdx,:);

            N = numel(xEstHouse);

            score = electionInfo.Score;
            geographyType = electionInfo.GeographyType;

            % Deal with uncontested races
            noDemFlag = electionInfo.CandidateD == "NaN";
            noRepFlag = electionInfo.CandidateR == "NaN";
            xEstHouse(noDemFlag) = 0;
            xEstHouse(noRepFlag) = 1;
            xSimsHouse(noDemFlag, :) = 0;
            xSimsHouse(noRepFlag, :) = 1;

            % Get seats per run
            nSims = size(xSimsHouse, 2);
            houseSeats = score' * (xSimsHouse > 0.5);
            houseMean = mean(houseSeats);
            houseWin = sum(houseSeats >= 218) / nSims;
            houseLose = sum(houseSeats < 218) / nSims;

            % Get popular vote
            hDistrict2Popular = 0 * xEstHouse;
            idxDistrict = geographyType == "Congressional District";
            prevVoteTot = model.districtData.districtInfo.TotalVote;
            H = model.createMapping(geographyType, electionInfo.GeographyName);
            hDistrict2Popular(idxDistrict) = H * prevVoteTot;
            hDistrict2Popular = hDistrict2Popular / sum(hDistrict2Popular);
            meanPopVote = xEstHouse' * hDistrict2Popular;
            
            sigmaPopVote = sqrt(hDistrict2Popular' * covEstHouse * ...
                hDistrict2Popular);

            % Display the results

            disp("--------------------------------------------------------")
            disp("House Election:")
            disp("  Democrats: " + num2str(houseMean,"%.2f") + " Seats | Chance of Winning: " + num2str(houseWin*100, "%.2f") + "%");
            disp("  Republicans: " + num2str(435-houseMean,"%.2f") + " Seats | Chance of Winning: " + num2str(houseLose*100, "%.2f") + "%");
            disp("--------------------------------------------------------")

            disp("House Popular Vote:")
            popVoteChance = normcdf(meanPopVote, 0.5, sigmaPopVote);
            disp("  Democrats: " + num2str(meanPopVote*100,"%.2f") + "% | Chance of Winning: " + num2str(popVoteChance * 100,"%.2f") + "%");
            disp("  Republicans: " + num2str(100-meanPopVote*100,"%.2f") + "% | Chance of Winning: " + num2str(100-popVoteChance * 100,"%.2f") + "%");
            disp("--------------------------------------------------------")

            electionNames = electionInfo.ElectionName;
            filterFlag = model.filterFlag(houseIdx);
            for i = 1:N-1
                xEstCur = xEstHouse(i);
                dCand = electionInfo.CandidateD(i);
                rCand = electionInfo.CandidateR(i);
                electionName = electionNames(i);
                if filterFlag(i)
                    sigmaEst = sqrt(covEstHouse(i,i));
                    chance = normcdf(xEstCur, 0.5, sigmaEst);
                    disp(electionName + ":")
                    disp("  " + dCand + ": " + num2str(xEstCur*100,"%.2f") + "% | Chance of Winning: " + num2str(chance*100,"%.2f") + "%")
                    disp("  " + rCand + ": " + num2str((1-xEstCur)*100,"%.2f") + "% | Chance of Winning: " + num2str((1-chance)*100,"%.2f") + "%")
                    disp(" ")
                else
                    if xEstCur == 0
                        disp(electionName + ":")
                        disp("  " + rCand + ": " + num2str(100,"%.2f") + "% | Chance of Winning: " + num2str(100,"%.2f") + "%")
                        disp(" ")
                    else
                        disp(electionName + ":")
                        disp("  " + dCand + ": " + num2str(100,"%.2f") + "% | Chance of Winning: " + num2str(100,"%.2f") + "%")
                        disp(" ")
                    end
                end
            end
        end

        function printSenateResults(model, xSims, ecVotes)
            arguments
                model matlab.Core.Model
                xSims (:,:) double
                ecVotes (1,:) double
            end

            % Model Senate Results
            senateIdx = model.electionInfo.ElectionType == "Senate";
            electionInfo = model.electionInfo(senateIdx,:);
            xEstSenate = model.xEst(senateIdx,:);
            covEstSenate = model.covEst(senateIdx,senateIdx);
            xSimsSenate = xSims(senateIdx,:);

            N = numel(xEstSenate);
            score = electionInfo.Score;

            % Deal with uncontested races
            noDemFlag = electionInfo.CandidateD == "NaN";
            noRepFlag = electionInfo.CandidateR == "NaN";
            xEstSenate(noDemFlag) = 0;
            xEstSenate(noRepFlag) = 1;
            xSimsSenate(noDemFlag, :) = 0;
            xSimsSenate(noRepFlag, :) = 1;

            % Get seats per run
            nSims = size(xSimsSenate, 2);
            senateSeats = score' * (xSimsSenate > 0.5) + model.config.senate.seats.dem;
            senateMean = mean(senateSeats);
            senateWin = (sum(senateSeats > 50) + sum(senateSeats == 50 & ecVotes >= 270)) / nSims;
            senateLose = (sum(senateSeats < 50) + sum(senateSeats == 50 & ecVotes < 270)) / nSims;

            % Display the results

            disp("--------------------------------------------------------")
            disp("Senate:")
            disp("  Democrats: " + num2str(senateMean,"%.2f") + " Seats | Chance of Winning: " + num2str(senateWin*100, "%.2f") + "%");
            disp("  Republicans: " + num2str(100-senateMean,"%.2f") + " Seats | Chance of Winning: " + num2str(senateLose*100, "%.2f") + "%");
            disp("--------------------------------------------------------")

            electionNames = electionInfo.ElectionName;
            filterFlag = model.filterFlag(senateIdx);
            for i = 1:N-1
                xEstCur = xEstSenate(i);
                dCand = electionInfo.CandidateD(i);
                rCand = electionInfo.CandidateR(i);
                electionName = electionNames(i);
                if filterFlag(i)
                    sigmaEst = sqrt(covEstSenate(i,i));
                    chance = normcdf(xEstCur, 0.5, sigmaEst);
                    disp(electionName + ":")
                    disp("  " + dCand + ": " + num2str(xEstCur*100,"%.2f") + "% | Chance of Winning: " + num2str(chance*100,"%.2f") + "%")
                    disp("  " + rCand + ": " + num2str((1-xEstCur)*100,"%.2f") + "% | Chance of Winning: " + num2str((1-chance)*100,"%.2f") + "%")
                    disp(" ")
                else
                    if xEstCur == 0
                        disp(electionName + ":")
                        disp("  " + rCand + ": " + num2str(100,"%.2f") + "% | Chance of Winning: " + num2str(100,"%.2f") + "%")
                        disp(" ")
                    else
                        disp(electionName + ":")
                        disp("  " + dCand + ": " + num2str(100,"%.2f") + "% | Chance of Winning: " + num2str(100,"%.2f") + "%")
                        disp(" ")
                    end
                end
            end
        end

        function printGubernatorialResults(model)
            arguments
                model matlab.Core.Model
            end

            % Model Gubernatorial Results
            govIdx = model.electionInfo.ElectionType == "Gubernatorial";
            electionInfo = model.electionInfo(govIdx,:);
            xEstGov = model.xEst(govIdx,:);
            covEstGov = model.covEst(govIdx,govIdx);

            N = numel(xEstGov);

            % Deal with uncontested races
            noDemFlag = electionInfo.CandidateD == "NaN";
            noRepFlag = electionInfo.CandidateR == "NaN";
            xEstGov(noDemFlag) = 0;
            xEstGov(noRepFlag) = 1;

            % Display the results

            disp("--------------------------------------------------------")
            disp("Gubernatorial Estimates:")
            disp("--------------------------------------------------------")

            electionNames = electionInfo.ElectionName;
            filterFlag = model.filterFlag(govIdx);
            for i = 1:N-1
                xEstCur = xEstGov(i);
                dCand = electionInfo.CandidateD(i);
                rCand = electionInfo.CandidateR(i);
                electionName = electionNames(i);
                if filterFlag(i)
                    sigmaEst = sqrt(covEstGov(i,i));
                    chance = normcdf(xEstCur, 0.5, sigmaEst);
                    disp(electionName + ":")
                    disp("  " + dCand + ": " + num2str(xEstCur*100,"%.2f") + "% | Chance of Winning: " + num2str(chance*100,"%.2f") + "%")
                    disp("  " + rCand + ": " + num2str((1-xEstCur)*100,"%.2f") + "% | Chance of Winning: " + num2str((1-chance)*100,"%.2f") + "%")
                    disp(" ")
                else
                    if xEstCur == 0
                        disp(electionName + ":")
                        disp("  " + rCand + ": " + num2str(100,"%.2f") + "% | Chance of Winning: " + num2str(100,"%.2f") + "%")
                        disp(" ")
                    else
                        disp(electionName + ":")
                        disp("  " + dCand + ": " + num2str(100,"%.2f") + "% | Chance of Winning: " + num2str(100,"%.2f") + "%")
                        disp(" ")
                    end
                end
            end
        end
    end
end