function Sadj = correct_sqrt_error(Sraw);
% This code post-corrects salinity measurements to account for a
% discontinuity in onboard salinity computation

% VARIABLES:
% Sraw - measured salinity
% Sadj - adjusted salinity

% Adjusts salinity
Sadj = Sraw;

% Set correction flag
correction = false;

% Loops throug hthe time series
for ii = 1:length(Sraw)

    if correction == false
        if Sraw(ii) < 35.000
            % Update correction flag
            correction = true;
        end
    else
        if Sraw(ii) > 35.002
            % Update correction flag
            correction = false;
        end
    end
    if correction == true
        % computes the error
        error = (3.559e-10)*exp(0.4403*Sraw(ii));
        % Adjusts the raw salinity
        Sadj(ii) = Sraw(ii) - error;
    end
end