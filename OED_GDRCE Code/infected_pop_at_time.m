function infpop = infected_pop_at_time(x,state,time)
% send times you want to get no. infected
infpop(1:length(x))=zeros;
for i=1:length(x)
    inf_at_timei=(state(time<=x(i)));
    infpop(i)=inf_at_timei(end);
end