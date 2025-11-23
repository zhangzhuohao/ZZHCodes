function fp_rand = gen_rand_fp(fp_target, fp_prob)

fp_pd = fp_prob ./ sum(fp_prob);
fp_cd = cumsum(fp_pd);

rand_i = rand();

fp_ind = find(rand_i<=fp_cd, 1, "first");
if isempty(fp_ind)
    fp_ind = length(fp_target);
end

fp_rand = fp_target(fp_ind);

end