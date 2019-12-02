function gamma=gaussian_model(h, p)
% gamma=gaussian_model(h, p)
% Compute a gaussian model, defined as
%         
% gamma= c{1 - exp(-h^2/a^2)}
% where
%	h is the lag
%	p(1) is the range
%	p(2) is the sill
%
% 

gamma =  p(2).*(1 - exp(-(h./p(1)).^2));



