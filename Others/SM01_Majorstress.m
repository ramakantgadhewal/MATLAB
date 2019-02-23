function [sigma, sigma_r] = SM01_Majorstress(S_prm)
    sigma(1,1) = 0.5 * (S_prm(1) + S_prm(2)) + 0.5 * sqrt((S_prm(1) - S_prm(2))^2 ...
        + 4 * S_prm(3)^2);
    sigma(2,1) = 0.5 * (S_prm(1) + S_prm(2)) - 0.5 * sqrt((S_prm(1) - S_prm(2))^2 ...
        + 4 * S_prm(3)^2);
    sigma_temp = [abs(sigma(1,1)), abs(sigma(2,1)), abs(sigma(1,1)-sigma(2,1))];
    sigma_r(1,1) = max(sigma_temp);
    sigma_r(2,1) = sqrt(sigma(1,1)^2 - sigma(1,1)*sigma(2,1) + sigma(2,1)^2);
end