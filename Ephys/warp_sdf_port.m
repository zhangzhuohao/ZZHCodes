function [sdf_warp_port, t_warped_port, t_points_port] = warp_sdf_port(sdf_warp, t_warped, t_points)

dt = unique(diff(t_warped{1}));
t_temp = cellfun(@(t_p, t_w) [t_w(1), t_p, t_w(end)+dt], t_points, t_warped, 'UniformOutput', false);

t_temp_m = cellfun(@(x1, x2) (x1 + x2) / 2, t_temp(:,1), t_temp(:,2), 'UniformOutput', false);
t_warp_m = cellfun(@(x) x(1):dt:x(end)-dt, t_temp_m, 'UniformOutput', false);

sdf_warp_port = cell(size(sdf_warp));
for i = 1:size(sdf_warp, 1)
    for j = 1:size(sdf_warp, 2)
        sdf_warp_ij = cell(1,length(t_temp_m{i})-1);
        for k = 1:length(t_temp_m{i})-1
            % template
            t_target_k = t_warp_m{i}(t_warp_m{i}>=t_temp_m{i}(k) & t_warp_m{i}<t_temp_m{i}(k+1));

            % origin
            t_towarp_k = t_warped{i,j}(t_warped{i,j}>=t_temp{i,j}(k) & t_warped{i,j}<t_temp{i,j}(k+1));
            sdf_towarp_k = sdf_warp{i,j}(:,t_warped{i,j}>=t_temp{i,j}(k) & t_warped{i,j}<t_temp{i,j}(k+1));

            % warp
            sdf_warp_ij{k} = warp_sdf_pop(t_towarp_k, sdf_towarp_k, t_target_k);
        end
        sdf_warp_port{i,j} = cat(2, sdf_warp_ij{:});
    end
end

t_warped_port = t_warp_m;
t_points_port = cellfun(@(t_p) t_p(2:end-1), t_temp_m, 'UniformOutput', false);

end

