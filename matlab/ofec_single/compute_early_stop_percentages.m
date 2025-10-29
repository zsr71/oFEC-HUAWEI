function tile_pct = compute_early_stop_percentages(tile_stats)
%COMPUTE_EARLY_STOP_PERCENTAGES Convert counters to percentages.
tile_pct = zeros(1, numel(tile_stats));
for i = 1:numel(tile_stats)
    total = tile_stats(i).total;
    if total > 0
        tile_pct(i) = double(tile_stats(i).triggered) / double(total) * 100.0;
    else
        tile_pct(i) = 0.0;
    end
end
end
