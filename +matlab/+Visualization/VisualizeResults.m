classdef VisualizeResults
    methods (Static)
        function printResults(model, xSims)
            arguments
                model matlab.Core.Model
                xSims (:,:) double
            end

            matlab.Visualization.VisualizeResults.printPresidentialResults(model, xSims);
        end

        function printPresidentialResults(model, xSims)
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
            N = numel(xEstPres(:,end));
            hX2Y = zeros(N, N - 1);
            hX2Y(~idxNat, :) = eye(N-1);
            hX2Y(idxNat, :) = ones(1, N-1);

            xEst = hX2Y' * xEstPres(:,end);
            xSimsPres = hX2Y' * xSimsPres;
            covEst = hX2Y' * covEstPres(:,:,end) * hX2Y;
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
            if electionInfo.IncumbencyFlag(idxNat) == 1
                dCand = dCand + "* (D)";
                rCand = rCand + " (R)";
            elseif electionInfo.IncumbencyFlag(idxNat) == -1
                dCand = dCand + " (D)";
                rCand = rCand + "* (R)";
            else
                dCand = dCand + " (D)";
                rCand = rCand + " (R)";
            end

            disp("Presidential Election:")
            disp("  " + dCand + ": " + num2str(ecMean,"%.2f") + " Electoral Votes | Chance of Winning: " + num2str(presECWin*100, "%.2f") + "%");
            disp("  " + rCand + ": " + num2str(538-ecMean,"%.2f") + " Electoral Votes | Chance of Winning: " + num2str(presECLose*100, "%.2f") + "%");
            disp("  Chance of Tie: " + num2str(presECTie*100,"%.2f") + "%");
            disp("")

            disp("Popular Vote:")
            popVoteChance = normcdf(meanPopVote, 0.5, sigmaPopVote);
            disp("  " + dCand + ": " + num2str(meanPopVote*100,"%.2f") + "% | Chance of Winning: " + num2str(popVoteChance * 100,"%.2f") + "%");
            disp("  " + rCand + ": " + num2str(100-meanPopVote*100,"%.2f") + "% | Chance of Winning: " + num2str(100-popVoteChance * 100,"%.2f") + "%");

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
                disp("")
                disp(geoName + " (" + num2str(score(i),"%.0f") + " Electoral Votes):")
                disp("  " + dCand + ": " + num2str(xEstCur*100,"%.2f") + "% | Chance of Winning: " + num2str(chance*100,"%.2f") + "%")
                disp("  " + rCand + ": " + num2str((1-xEstCur)*100,"%.2f") + "% | Chance of Winning: " + num2str((1-chance)*100,"%.2f") + "%")
            end
        end
    end
end