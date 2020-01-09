%% Poisson distribution
mu = 4;
poiss_RNG = poissrnd(mu, [1000, 1]);
csvwrite('poiss_RNG.csv', poiss_RNG)

%% Power distribution
upper = 6;
lower = 2;
shape = 1;    % Shape parameter
scale = 4;   % Scale parameter
pow_RNG = gamrnd(shape, scale, [1000, 1]);  % Generate 1000 random numbers with gamma distribution
% pow_RNG = exp(pow_RNG);

% range = max(pow_RNG) - min(pow_RNG);
% normPow_RNG = (pow_RNG - min(pow_RNG))/range;
% 
% range2 = upper - lower;
% normPow_RNG = (normPow_RNG*range2) + lower;
histogram(pow_RNG)
csvwrite('pow_RNG.csv', pow_RNG)  % Write to file

%% Uniform distribution
lower = 2;
upper = 6;
uniform_RNG = unifrnd(lower, upper, [1000 1]);

csvwrite('Uniform_RNG.csv', uniform_RNG)  % Write to file

%% Analysis as in Hamill

