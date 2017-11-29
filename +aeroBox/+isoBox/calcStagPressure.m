function P0 = calcStagPressure(varargin)
%

np = inputParser();
np.addParameter('mach', @isnumeric);
np.addParameter('gamma', @isnumeric);
np.addParameter('Ps', @isnumeric);
np.parse(varargin{:});

M = np.Results.mach;
gamma = np.Results.gamma;
P = np.Results.Ps;

P0 = P * ((1 + (((gamma - 1) / 2) * M^2))^(gamma / (gamma - 1)));


end

