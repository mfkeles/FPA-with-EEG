function score_arr = find_scores(fiber_time,scores)
%Finds scores for the fiber data based on the time from start in the score table

positions = interp1(scores.TimeFromStart, 1:numel(scores.TimeFromStart), fiber_time);

score_arr = scores.Qiang_0_Numeric(floor(positions));

if find(score_arr == 129)
    score_arr(score_arr==129) = 1;
end

