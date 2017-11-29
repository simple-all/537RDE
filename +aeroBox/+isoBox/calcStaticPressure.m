function P = calcStaticPressure(varargin)

np = inputParser();
np.addParameter('mach', @isnumeric);
np.addParameter('gamma', @isnumeric);
np.addParameter('Pt', @isnumeric);
np.parse(varargin{:});

M = np.Results.mach;
gamma = np.Results.gamma;
P_0 = np.Results.Pt;

P = P_0 / ((1 + (((gamma - 1) / 2) * M^2))^(gamma / (gamma - 1)));


end
