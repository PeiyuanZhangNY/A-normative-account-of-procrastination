function [patternLabels, patternCounts, patternNames, patternMatrices] = classifyTimeCoursePatterns(WithoptActSeqMatrix3)
% classifyTimeCoursePatterns  Classify each 10-day work-amount time course
% into one of six qualitative shape patterns.
%
%   [patternLabels, patternCounts, patternNames, patternMatrices] = ...
%       classifyTimeCoursePatterns(WithoptActSeqMatrix3)
%
% INPUT
%   WithoptActSeqMatrix3 : N x T matrix (e.g. 10000 x 10). Each row is a
%                          time course of work amount over T days.
%
% OUTPUT
%   patternLabels   : N x 1 integer label (1..6) for each row.
%   patternCounts   : 6 x 1 count of rows in each pattern.
%   patternNames    : 1 x 6 cell array of the pattern names.
%   patternMatrices : 1 x 6 cell array. patternMatrices{k} is an
%                     (n_k x T) matrix containing every time course
%                     assigned to pattern k, in their original row order.
%                     Empty patterns return a 0 x T matrix.
%
% Tip: to also recover which original rows went into each pattern, use
% patternLabels, e.g. rowsRampDown = find(patternLabels == 3).
%
% PATTERN CODES
%   1 = Steady          (roughly flat)
%   2 = Ramping up      (monotone-ish increasing)
%   3 = Going down      (monotone-ish decreasing)
%   4 = Up then down    (inverted-U / single interior peak)
%   5 = Down then up    (U-shape / single interior trough)
%   6 = Fluctuating     (multiple direction changes, none of the above)

    [N, T] = size(WithoptActSeqMatrix3);

    patternNames = {'Steady','RampUp','RampDown','UpThenDown','DownThenUp','Fluctuating'};
    patternLabels = zeros(N,1);

    % --- Tunable thresholds -------------------------------------------------
    % "Steady" tolerance: total range relative to the mean level of the row.
    steadyRangeFrac = 0.10;   % range < 10% of mean level => steady
    % Minimum slope (per step) to count a segment as a real trend, expressed
    % as a fraction of the row's range. Small wiggles are treated as noise.
    noiseFrac       = 0.05;
    % ------------------------------------------------------------------------

    for i = 1:N
        x = WithoptActSeqMatrix3(i,:);
        rng = max(x) - min(x);
        mu  = mean(abs(x));

        % ----- 1. Steady: essentially flat -----
        if rng < steadyRangeFrac * max(mu, eps)
            patternLabels(i) = 1;
            continue;
        end

        % Noise floor in absolute units for this row
        noiseThr = noiseFrac * rng;

        % First differences, then sign with a dead zone for noise
        d = diff(x);
        s = zeros(1, T-1);
        s(d >  noiseThr) =  1;
        s(d < -noiseThr) = -1;

        % Drop zeros (flat sub-segments) and merge consecutive equal signs
        s = s(s ~= 0);
        if isempty(s)
            % Range exceeded steady threshold but no segment cleared the
            % noise floor: treat as steady.
            patternLabels(i) = 1;
            continue;
        end
        s = s([true, diff(s) ~= 0]);   % run-length encode sign changes

        % Number of monotone segments after merging
        nSeg = numel(s);

        if nSeg == 1
            if s == 1
                patternLabels(i) = 2;   % ramp up
            else
                patternLabels(i) = 3;   % ramp down
            end
        elseif nSeg == 2
            if isequal(s, [1 -1])
                patternLabels(i) = 4;   % up then down
            elseif isequal(s, [-1 1])
                patternLabels(i) = 5;   % down then up
            else
                patternLabels(i) = 6;
            end
        else
            patternLabels(i) = 6;       % >=3 turns => fluctuating
        end
    end

    patternCounts = accumarray(patternLabels, 1, [6 1]);

    % --- Build one matrix of time courses per pattern -----------------------
    % patternMatrices{k} contains every row of the input whose label is k.
    patternMatrices = cell(1, 6);
    for k = 1:6
        patternMatrices{k} = WithoptActSeqMatrix3(patternLabels == k, :);
    end

    % --- Report -------------------------------------------------------------
    fprintf('\nPattern classification of %d time courses (T = %d days):\n', N, T);
    fprintf('%-14s %8s %8s\n', 'Pattern', 'Count', 'Percent');
    for k = 1:6
        fprintf('%-14s %8d %7.1f%%\n', patternNames{k}, ...
                patternCounts(k), 100*patternCounts(k)/N);
    end
    fprintf('Distinct pattern kinds present: %d\n\n', nnz(patternCounts));
end
