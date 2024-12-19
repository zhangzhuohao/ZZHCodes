function r = shiftWaldRanNum(alpha,theta,gamma,n,m)
 
r = zeros(n,m);
for i = 1:n

    % Generate a random value from a normal distribution with mean 0 
    % and standard deviation equal to 1
    y = normrnd(0,1);

    % Square the value
    y = y^2;

    % Now use the relation
    x = alpha+(alpha^2*y)/(2*gamma)-alpha/(2*gamma) * ...
        sqrt((4*alpha*gamma*y)+alpha^2*y^2);

    % Generate another random variable, this time sampled from a uniform
    % distribution between 0 and 1
    z = rand(1);
    
    % Finally ...
    if (z < alpha/(alpha+x))
        r(i) = theta + x; 
    else
        r(i) = theta + (alpha^2/x); 
    end
    
end