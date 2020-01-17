function generateDistributions(distributionType, distributionParameters)
    switch distributionType
        case 'poisson'
            % ----------------------Poisson distribution----------------------
            poiss_RNG = poissrnd(distributionParameters, [1000, 1]);
%             csvwrite('poiss_RNG.csv', poiss_RNG)
            csvwrite('randomDistribution.csv', poiss_RNG)
        
        case 'power'
            % ----------------------Power distribution----------------------
            lower = distributionParameters(1);
            upper = distributionParameters(2);
            mu = distributionParameters(3);
            sigma = distributionParameters(4);
            
            shape = mu^2/sigma^2;
            scale = sigma^2/mu;
            
            pow_RNG = gamrnd(shape, scale, [1000, 1]);  % Generate 1000 random numbers with gamma distribution
            pow_RNG = exp(pow_RNG);

            pow_RNG(pow_RNG > upper) = upper;
            pow_RNG(pow_RNG < lower) = lower;
%             csvwrite('pow_RNG.csv', pow_RNG)  % Write to file
            csvwrite('randomDistribution.csv', pow_RNG)
        
        case 'uniform'
            % ----------------------Uniform distribution----------------------
            lower = distributionParameters(1);
            upper = distributionParameters(2);
            uniform_RNG = unifrnd(lower, upper, [1000 1]);
            csvwrite('randomDistribution.csv', uniform_RNG)
%             csvwrite('Uniform_RNG.csv', uniform_RNG)  % Write to file
    end
end


