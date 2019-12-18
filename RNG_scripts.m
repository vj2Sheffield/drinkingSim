%% Poisson distribution
mu = 22;
poiss_RNG = poissrnd(mu, [1000, 1]);
csvwrite('poiss_RNG.csv', poiss_RNG)

%% Power distribution

clear, clc
shape = 400;    % Shape parameter
scale = 1/20;   % Scale parameter
pow_RNG = gamrnd(shape, scale, [1000, 1]);  % Generate 1000 random numbers with gamma distribution

csvwrite('pow_RNG.csv', pow_RNG)  % Write to file




