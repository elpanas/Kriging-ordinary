function gamma=spherical_model(h, p)
% gamma=spherical_model(h, p)
% Compute a spherical model, defined as
%         |0 if h=0
% gamma= < c*(3h/(2a) - 1/2*(h/a)^3)
%         |c_0 if h>a
% where
%	h is the lag
%	p(1) is the range
%	p(2) is the  sill
% 

gamma =  p(2).*((3*h)./(2*p(1)) - (0.5 .* (h./p(1)).^3));

idx = find (h > p(1));
gamma(idx) = p(2);


