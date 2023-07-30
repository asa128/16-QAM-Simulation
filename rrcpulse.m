function [h,t] = rrcpulse(k_T, T_sym, s, r)
%SINCPULSE Generates a root raised cosine rolloff pulse with 2*k_T symbol periods, symbol period T_sym, s samples per symbol, and rolloff factor r.
%   Returns a vector of pulse samples and a vector of corresponding times

D = 1 / T_sym; % symbol rate

t_step = T_sym/s; % time step

t = -k_T*T_sym : t_step : ((k_T)*T_sym - t_step); % array of times

% RRC pulse function
h = ( sin(pi.*D.*t.*(1-r)) + 4.*D.*r.*t.*cos(pi.*D.*t.*(1+r)) ) ./ ( pi.*D.*t.*(1-(4.*D.*r.*t).^2) );

% address special cases where the above function is undefined
k = T_sym/(4*r);
h( abs(t-k) < 0.1*t_step ) = (r./sqrt(2)).*( (1+2./pi).*sin(pi./(4.*r)) + (1-2./pi).*cos(pi./(4.*r)) );
h( abs(t+k) < 0.1*t_step ) = (r./sqrt(2)).*( (1+2./pi).*sin(pi./(4.*r)) + (1-2./pi).*cos(pi./(4.*r)) );
h( abs(t) < 0.1*t_step ) = 1 - r + 4*r/pi;

end