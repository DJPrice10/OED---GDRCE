function new_design = discrete_trunc_sample(T, minmaxT, sd_T, stepsize, m)
% % Sample new design points corresponding to constraints, and rounded to a grid 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

new_design = -1*ones(m,length(T));
for i=1:m
    while( any( any(new_design(i,:)<min(minmaxT)) || any(new_design(i,:)>max(minmaxT)) )  )
        new_design(i,:) = round( mvnrnd(T, sd_T, 1) /stepsize) * stepsize;
    end
end
