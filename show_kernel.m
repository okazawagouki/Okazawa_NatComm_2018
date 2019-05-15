function show_kernel(sim, alignment_type, smoothing_wid, xrange, yrange)


kernel_len = 1000;
weight_len = xrange(2);

smoothen = @(x, val) conv(x, fspecial('average',[1 val]),'same')./conv(ones(size(x)), fspecial('average',[1 val]),'same');


bound_height = mean(sim.param.bound_height(1:round(sim.medianRT), 2)) - sim.param.B0;

scaling_factor = 1 / bound_height * 2;

switch alignment_type
    case 'stim'
        kernel = sim.kernel.stim{1} - sim.kernel.stim{2};
    case 'resp'
        kernel = sim.kernel.resp{1} - sim.kernel.resp{2};
end
kernel = smoothen(kernel, smoothing_wid);

kernel = kernel / scaling_factor;
weight = sim.param.weight(1:weight_len);


hold on;

switch alignment_type
    case 'stim'
        plot(0:(length(weight)-1), weight, '--', 'color', [.3 .3 .3], 'linew', 1);
        plot(0:kernel_len-1, kernel(1:kernel_len), 'linew', 2, 'color', 'k');
        set(gca, 'xlim', [0 1200], 'xtick', 0:250:1000, 'xticklabel', [], 'ylim', yrange([1 end]), 'ytick', yrange, ...
             'box', 'off', 'tickdir', 'out', 'ticklen', [.03 .03]);
        axis square;
        ylabel('Weight', 'fontsize', 10);
        set(gca, 'xticklabel', {'0', '', '500', '', '1000'});
        xlabel({'Time from', 'stimulus onset (ms)'}, 'fontsize', 10);
    case 'resp'
        plot((-kernel_len+1):1:0, kernel(end-kernel_len+1:end), 'linew', 2, 'color', 'k');
        set(gca, 'xlim', [-1200 -20], 'xtick', -1000:250:0, 'xticklabel', [], 'ylim', yrange([1 end]), 'ytick', yrange,...
             'YAxisLocation', 'right', 'box', 'off', 'tickdir', 'out', 'ticklen', [.03 .03]);
        axis square;
        set(gca, 'xticklabel',{'-1000', '', '-500', '', '0'});
        xlabel({'Time from', 'response onset (ms)'}, 'fontsize', 10);
end

