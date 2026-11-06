% ============================================================
% MAIN CWD SIMULATION
% ============================================================
%
% Spatial, sex-structured CWD SIH model based on:
%
% Garlick et al. (2014)
% "Homogenization, sex, and differential motility"
%
% HOMOGENIZED FORMULATION
%
% Model variables:
%
%   cIF = homogenized infected-female variable
%   cIM = homogenized infected-male variable
%   H   = environmental hazard
%
% Physical infected densities:
%
%   IF = cIF ./ MuFemale
%   IM = cIM ./ MuMale
%
% Susceptible steady-state approximation:
%
%   SF = CF ./ MuFemale
%   SM = CM ./ MuMale
%
% The constants CF and CM are selected so that:
%
%   integral(SF) = NF
%   integral(SM) = NM
%
%
% SPATIAL RESOLUTION
%
%   LANDFIRE              ~30 m
%   Model cells           600 m
%   Homogenization block  3 km
%
%
% SEASONALITY
%
%   Summer    = 153 days
%   Breeding  = 61 days
%   Winter    = 151 days
%
%
% INITIAL INFECTION
%
% The paper places the initial CWD infection in the
% La Sal Mountains.
%
% The infection is placed directly on the user's LANDFIRE
% raster using the approximate location marked on the TIF:
%
%   approximately 90.5% across the raster
%   approximately 24.3% down from the top
%
% This avoids requiring geographic-coordinate lookup.
%
% Initial infected population:
%
%   47 infected females
%   32 infected males
%
% The infection is distributed as a small Gaussian hotspot
% around the approximate La Sal Mountains location.
%
% The infected population and environmental hazard are NOT
% reset when seasons change.
%
% ============================================================


clear;
clc;


%% ============================================================
% 1. INPUT LANDFIRE RASTER
% ============================================================

tifFile = fullfile( ...
    'data', ...
    'LF2025_EVT_cropped5percent.tif');


if ~isfile(tifFile)

    error( ...
        'LANDFIRE file not found: %s', ...
        tifFile);

end


%% ============================================================
% 2. MODEL PARAMETERS
% ============================================================

params = struct();


%% ------------------------------------------------------------
% Population
% ------------------------------------------------------------

params.NF = 22430;
params.NM = 4070;

params.Ntotal = ...
    params.NF + params.NM;


%% ------------------------------------------------------------
% Spatial resolution
% ------------------------------------------------------------

params.dx = 0.6;
params.dy = 0.6;

params.cellArea = ...
    params.dx * params.dy;


%% ============================================================
% 3. SIMULATION TIME
% ============================================================

% Ten years

params.T = 3650;

% Garlick et al. time step

params.dt = 5.0;

% Save every 5 days

params.outputDt = 5.0;


%% ============================================================
% 4. SEASONALITY
% ============================================================

params.seasonNames = { ...
    'summer', ...
    'breeding', ...
    'winter' };


params.seasonLength = [ ...
    153, ...
     61, ...
    151 ];


params.yearLength = ...
    sum(params.seasonLength);


%% ============================================================
% 5. DIRECT TRANSMISSION
% ============================================================

params.betaFF = 0.000625;
params.betaMF = 0.000625;

params.betaFM = 0.0019;
params.betaMM = 0.0019;


% ------------------------------------------------------------
% Summer
% ------------------------------------------------------------

params.betaFF_summer = ...
    params.betaFF;

params.betaMF_summer = ...
    params.betaMF;

params.betaFM_summer = ...
    params.betaFM;

params.betaMM_summer = ...
    params.betaMM;


% ------------------------------------------------------------
% Breeding
% ------------------------------------------------------------

params.betaFF_breeding = ...
    params.betaFF;

params.betaMF_breeding = ...
    params.betaMF;

params.betaFM_breeding = ...
    params.betaFM;

params.betaMM_breeding = ...
    params.betaMM;


% ------------------------------------------------------------
% Winter
% ------------------------------------------------------------

params.betaFF_winter = ...
    params.betaFF;

params.betaMF_winter = ...
    params.betaMF;

params.betaFM_winter = ...
    params.betaFM;

params.betaMM_winter = ...
    params.betaMM;


%% ============================================================
% 6. ENVIRONMENTAL TRANSMISSION
% ============================================================

params.gammaF = 0.0022;
params.gammaM = 0.0022;


%% ============================================================
% 7. ENVIRONMENTAL CONTAMINATION
% ============================================================

params.alphaF = 0.00030411;
params.alphaM = 0.00030411;

params.alphaDF = 0.1;
params.alphaDM = 0.1;


%% ============================================================
% 8. CWD MORTALITY
% ============================================================

params.omega = 0.0013;


%% ============================================================
% 9. ENVIRONMENTAL DECAY
% ============================================================

params.delta = 0.007;


%% ============================================================
% 10. CULLING
% ============================================================

params.phi = 0;


%% ============================================================
% 11. BUILD SPATIAL BACKGROUND
% ============================================================

[Landscape, ...
 ParSets, ...
 Mu, ...
 mubar, ...
 SS0A, ...
 SF_unused, ...
 SM_unused] = ...
    buildBackground( ...
        params, ...
        tifFile, ...
        params.NF, ...
        params.NM, ...
        params.cellArea);


%% ============================================================
% 12. CHECK REQUIRED SEASONAL FIELDS
% ============================================================

requiredFields = { ...
    'MuFemaleSummer', ...
    'MuFemaleBreeding', ...
    'MuFemaleWinter', ...
    'MuMaleSummer', ...
    'MuMaleBreeding', ...
    'MuMaleWinter' };


for k = 1:length(requiredFields)

    if ~isfield(ParSets,requiredFields{k})

        error( ...
            ['buildBackground did not return ' ...
             'the required seasonal field: %s'], ...
            requiredFields{k});

    end

end


%% ============================================================
% 13. STORE SEASONAL MOTILITY FIELDS
% ============================================================

params.MuFemaleSummer = ...
    ParSets.MuFemaleSummer;

params.MuFemaleBreeding = ...
    ParSets.MuFemaleBreeding;

params.MuFemaleWinter = ...
    ParSets.MuFemaleWinter;


params.MuMaleSummer = ...
    ParSets.MuMaleSummer;

params.MuMaleBreeding = ...
    ParSets.MuMaleBreeding;

params.MuMaleWinter = ...
    ParSets.MuMaleWinter;


%% ============================================================
% 14. INITIAL SEASON
% ============================================================

params.currentSeason = ...
    'summer';


params.MuFemale = ...
    params.MuFemaleSummer;

params.MuMale = ...
    params.MuMaleSummer;


%% ============================================================
% 15. CHECK MOTILITY
% ============================================================

if any(params.MuFemale(:) <= 0)

    error( ...
        'Summer female motility contains zero or negative values.');

end


if any(params.MuMale(:) <= 0)

    error( ...
        'Summer male motility contains zero or negative values.');

end


%% ============================================================
% 16. INITIAL SUSCEPTIBLE POPULATION
% ============================================================

[SF0,SM0,CF,CM] = ...
    initializeGarlickSusceptibles( ...
        params.MuFemale, ...
        params.MuMale, ...
        params.NF, ...
        params.NM, ...
        params.cellArea);


%% ============================================================
% 17. INITIAL POPULATION CHECK
% ============================================================

SF_initial = ...
    sum(SF0(:)) * params.cellArea;

SM_initial = ...
    sum(SM0(:)) * params.cellArea;


fprintf('\n');
fprintf('============================================\n');
fprintf(' GARLICK STEADY-STATE POPULATION\n');
fprintf('============================================\n');

fprintf('\nTarget population:\n');

fprintf(' Females: %.2f\n', ...
    params.NF);

fprintf(' Males:   %.2f\n', ...
    params.NM);

fprintf(' Total:   %.2f\n', ...
    params.Ntotal);


fprintf('\nCalculated summer population:\n');

fprintf(' Females: %.2f\n', ...
    SF_initial);

fprintf(' Males:   %.2f\n', ...
    SM_initial);

fprintf(' Total:   %.2f\n', ...
    SF_initial + SM_initial);


%% ============================================================
% 18. INITIAL INFECTED POPULATION
% ============================================================
%
% Garlick et al. simulations:
%
%   47 infected females
%   32 infected males
%
% ------------------------------------------------------------

IF_target = 47;
IM_target = 32;


%% ============================================================
% 19. APPROXIMATE LA SAL MOUNTAINS LOCATION
% ============================================================
%
% The paper identifies the La Sal Mountains as the location
% where the initial CWD infection is introduced.
%
% Instead of using geographic coordinates, the infection is
% placed directly using the location marked on the user's
% LANDFIRE TIF.
%
% Approximate position on the displayed TIF:
%
%   X = 90.5% across from the left
%   Y = 24.3% down from the top
%
% These values are an approximate visual placement based on
% the provided raster image. They are not geographic
% coordinates taken from the paper.
%
% ------------------------------------------------------------

laSalXFraction = 0.905;
laSalYFraction = 0.243;


%% ============================================================
% 20. CONVERT RASTER POSITION TO MODEL-GRID CELL
% ============================================================

[L1,L2] = ...
    size(SF0);


% Column increases from left to right

laSalCol = ...
    round(laSalXFraction * L2);


% Row increases from top to bottom in the MATLAB raster

laSalRow = ...
    round(laSalYFraction * L1);


% Make sure the selected location is inside the model grid

laSalRow = ...
    max(1,min(L1,laSalRow));

laSalCol = ...
    max(1,min(L2,laSalCol));


%% ============================================================
% 21. DISPLAY INITIAL INFECTION LOCATION
% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf(' INITIAL CWD LOCATION\n');
fprintf('============================================\n');

fprintf('\nLocation used from LANDFIRE raster:\n');

fprintf(' Horizontal position: %.1f%%\n', ...
    100 * laSalXFraction);

fprintf(' Vertical position:   %.1f%%\n', ...
    100 * laSalYFraction);

fprintf('\nModel-grid location:\n');

fprintf(' Row:    %d / %d\n', ...
    laSalRow, ...
    L1);

fprintf(' Column: %d / %d\n', ...
    laSalCol, ...
    L2);


%% ============================================================
% 22. INITIAL INFECTION DISTRIBUTION
% ============================================================
%
% Rather than putting all infected deer into one model cell,
% use a small Gaussian hotspot centered on the approximate
% La Sal Mountains location.
%
% sigma = 2 model cells
%
% At 600 m per model cell, this corresponds to a relatively
% small initial infection area centered on the La Sal range.
%
% ------------------------------------------------------------

[xGrid,yGrid] = ...
    meshgrid(1:L2,1:L1);


sigma = 2;


infectionMask = exp( ...
    -((xGrid-laSalCol).^2 + ...
      (yGrid-laSalRow).^2) ./ ...
      (2*sigma^2));


%% ============================================================
% 23. INITIAL INFECTED FEMALES
% ============================================================

IF0 = ...
    SF0 .* infectionMask;


% Normalize to exactly 47 infected females

IF0 = ...
    IF0 .* ...
    (IF_target / ...
    (sum(IF0(:)) * params.cellArea));


%% ============================================================
% 24. INITIAL INFECTED MALES
% ============================================================

IM0 = ...
    SM0 .* infectionMask;


% Normalize to exactly 32 infected males

IM0 = ...
    IM0 .* ...
    (IM_target / ...
    (sum(IM0(:)) * params.cellArea));


%% ============================================================
% 25. INITIAL ENVIRONMENTAL HAZARD
% ============================================================

H0 = ...
    zeros(size(SF0));


%% ============================================================
% 26. INITIAL CONDITION CHECK
% ============================================================

initialIF = ...
    sum(IF0(:)) * params.cellArea;

initialIM = ...
    sum(IM0(:)) * params.cellArea;


fprintf('\n');
fprintf('============================================\n');
fprintf(' INITIAL CONDITIONS\n');
fprintf('============================================\n');

fprintf('\nInitial season:\n');

fprintf(' %s\n', ...
    params.currentSeason);


fprintf('\nInitial infected population:\n');

fprintf(' Female: %.2f\n', ...
    initialIF);

fprintf(' Male:   %.2f\n', ...
    initialIM);

fprintf(' Total:  %.2f\n', ...
    initialIF + initialIM);


%% ============================================================
% 27. VERIFY INITIAL INFECTION LOCATION
% ============================================================
%
% Display the model landscape with the initial infection
% center marked.
%
% This is only a diagnostic plot and does not affect the
% simulation.
%
% ------------------------------------------------------------

figure;

imagesc(Landscape);

axis image;

hold on;

plot( ...
    laSalCol, ...
    laSalRow, ...
    'ro', ...
    'MarkerSize',12, ...
    'LineWidth',2);


title('Initial CWD Location - Approximate La Sal Mountains');

xlabel('Model column');
ylabel('Model row');

colorbar;

hold off;


%% ============================================================
% 28. SEASONAL SETTINGS
% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf(' SEASONAL MODEL\n');
fprintf('============================================\n');

fprintf('\nSummer:\n');
fprintf(' 153 days\n');

fprintf('Breeding:\n');
fprintf(' 61 days\n');

fprintf('Winter:\n');
fprintf(' 151 days\n');

fprintf('Annual cycle:\n');
fprintf(' 365 days\n');


%% ============================================================
% 29. MOTILITY RANGE
% ============================================================

fprintf('\n');
fprintf('Seasonal motility ranges:\n');


fprintf('\nFemale:\n');

fprintf(' Summer:   %.6g to %.6g\n', ...
    min(params.MuFemaleSummer(:)), ...
    max(params.MuFemaleSummer(:)));

fprintf(' Breeding: %.6g to %.6g\n', ...
    min(params.MuFemaleBreeding(:)), ...
    max(params.MuFemaleBreeding(:)));

fprintf(' Winter:   %.6g to %.6g\n', ...
    min(params.MuFemaleWinter(:)), ...
    max(params.MuFemaleWinter(:)));


fprintf('\nMale:\n');

fprintf(' Summer:   %.6g to %.6g\n', ...
    min(params.MuMaleSummer(:)), ...
    max(params.MuMaleSummer(:)));

fprintf(' Breeding: %.6g to %.6g\n', ...
    min(params.MuMaleBreeding(:)), ...
    max(params.MuMaleBreeding(:)));

fprintf(' Winter:   %.6g to %.6g\n', ...
    min(params.MuMaleWinter(:)), ...
    max(params.MuMaleWinter(:)));


%% ============================================================
% 30. RUN HOMOGENIZED GARLICK MODEL
% ============================================================

results = ...
    CWDSimulationGarlick( ...
        params, ...
        Landscape, ...
        ParSets, ...
        Mu, ...
        mubar, ...
        SF0, ...
        SM0, ...
        IF0, ...
        IM0, ...
        H0, ...
        CF, ...
        CM);


%% ============================================================
% 31. SAVE OUTPUTS
% ============================================================

saveOutputs( ...
    results, ...
    Landscape, ...
    tifFile);


%% ============================================================
% 32. FINAL SUMMARY
% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf(' MAIN SCRIPT COMPLETE\n');
fprintf('============================================\n');

fprintf('\nFinal infected females: %.4f\n', ...
    results.totalIF(end));

fprintf('Final infected males:   %.4f\n', ...
    results.totalIM(end));

fprintf('Final infected total:   %.4f\n', ...
    results.totalI(end));

fprintf('\nOutputs saved successfully.\n');


%% ============================================================
% LOCAL FUNCTION:
% INITIALIZE GARLICK STEADY-STATE SUSCEPTIBLES
% ============================================================

function [SF,SM,CF,CM] = ...
    initializeGarlickSusceptibles( ...
        MuFemale, ...
        MuMale, ...
        NF, ...
        NM, ...
        cellArea)


invMuF = ...
    1 ./ MuFemale;

invMuM = ...
    1 ./ MuMale;


CF = ...
    NF / ...
    (sum(invMuF(:)) * cellArea);


CM = ...
    NM / ...
    (sum(invMuM(:)) * cellArea);


SF = ...
    CF .* invMuF;


SM = ...
    CM .* invMuM;

end