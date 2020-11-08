function [log_whittle] = whittle_like_replicate(xobs, Beta)

nsubj = size(xobs,1);
log_whittle = 0;
for j=1:nsubj
   log_whittle = log_whittle + whittle_like_single(xobs{j},Beta); 
end    
