function [Landscape,ParSets,Mu,mubar,SS0A,SF,SM] = ...
    buildBackground(params,tifFile,NF,NM,cellArea)

% ============================================================
% buildBackground
% ============================================================
%
% Builds the spatial background for the Garlick et al. (2014)
% homogenized CWD SIH model.
%
% SPATIAL SCALE
%
%   LANDFIRE raster:       30 m
%   Model grid:           600 m
%   Homogenization scale: 3 km
%
% Therefore:
%
%   600 / 30 = 20
%   3000 / 600 = 5
%
% Processing:
%
%   1. Read 30-m LANDFIRE EVT raster
%   2. Convert LANDFIRE EVT -> Garlick EVT
%   3. Aggregate converted EVT to 600-m cells
%   4. Assign sex-specific motility
%   5. Calculate local 3-km harmonic means
%   6. Construct Garlick steady-state susceptible populations
%
% IMPORTANT
%
% The LANDFIRE -> Garlick conversion is performed BEFORE
% spatial aggregation. This ensures that multiple LANDFIRE
% classes mapping to the same Garlick class are treated as
% the same habitat before determining the dominant 600-m class.
%
% ============================================================


%% ============================================================
% 1. CHECK INPUTS
% ============================================================

if nargin ~= 5

    error([ ...
        'buildBackground requires exactly 5 inputs:\n' ...
        'buildBackground(params,tifFile,NF,NM,cellArea)' ]);

end


if ~isfile(tifFile)

    error( ...
        'LANDFIRE file not found: %s', ...
        tifFile);

end


%% ============================================================
% 2. SPATIAL RESOLUTION
% ============================================================

landfireResolution = 0.03;       % km = 30 m
modelResolution = 0.60;          % km = 600 m
homogenizationResolution = 3.00; % km = 3000 m


aggregationFactor = ...
    round(modelResolution / landfireResolution);


homogenizationFactor = ...
    round(homogenizationResolution / modelResolution);


if aggregationFactor < 1

    error('Invalid LANDFIRE aggregation factor.');

end


if abs( ...
        aggregationFactor * landfireResolution ...
        - modelResolution) > 1e-10

    error( ...
        'Model resolution must be an integer multiple of LANDFIRE resolution.');

end


if homogenizationFactor < 1

    error('Invalid homogenization factor.');

end


if abs( ...
        homogenizationFactor * modelResolution ...
        - homogenizationResolution) > 1e-10

    error( ...
        'Homogenization resolution must be an integer multiple of model resolution.');

end


fprintf('\n');
fprintf('============================================\n');
fprintf(' BUILDING SPATIAL BACKGROUND\n');
fprintf('============================================\n');

fprintf('\nSpatial resolution:\n');

fprintf('  LANDFIRE:       %.3f km\n', ...
    landfireResolution);

fprintf('  Model:          %.3f km\n', ...
    modelResolution);

fprintf('  Homogenization: %.3f km\n', ...
    homogenizationResolution);

fprintf('  Aggregation:    %d x %d LANDFIRE pixels\n', ...
    aggregationFactor,aggregationFactor);

fprintf('  Homogenization: %d x %d model cells\n', ...
    homogenizationFactor,homogenizationFactor);


%% ============================================================
% 3. READ LANDFIRE RASTER
% ============================================================

fprintf('\nLoading LANDFIRE raster...\n');


EVT = imread(tifFile);


if ndims(EVT) ~= 2

    error( ...
        'LANDFIRE raster must be a single-band raster.');

end


EVT = double(EVT);


[nRows,nCols] = size(EVT);


fprintf('  Original raster size: %d x %d\n', ...
    nRows,nCols);


%% ============================================================
% 4. COMPLETE 2025 LANDFIRE -> GARLICK CROSSWALK
% ============================================================
%
% First column:
%       2025 LANDFIRE EVT
%
% Second column:
%       Garlick EVT habitat class
%
% The additional classes at the bottom of this table are the
% classes detected in the user's actual raster that were not
% present in the previous crosswalk.
%
% ============================================================

crosswalk = [

    % ========================================================
    % OPEN / DEVELOPED
    % ========================================================

    7292    11

    7300    21
    7296    22

    7297    23
    7298    23
    7299    23

    7295    31

    7903    21
    7904    21


    % ========================================================
    % BARREN / AGRICULTURE
    % ========================================================

    9001    31
    9016    31
    9018    31
    9009    31

    7754    81
    7755    82

    7963    82
    7965    82
    7966    81
    7967    81
    7968    82


    % ========================================================
    % FOREST / WOODLAND
    % ========================================================

    7011    2011

    7016    2016
    7019    2016
    7059    2016
    7102    2016

    7051    2051

    7045    2052
    7047    2052
    7052    2052

    7054    2054
    7055    2055

    7049    2057
    7061    2061


    % ========================================================
    % SHRUBLANDS
    % ========================================================

    7062    2062
    7250    2062

    7064    2064
    7124    2064

    7066    2066

    7080    2080

    7075    2081
    7081    2081
    7088    2081

    7086    2086

    7093    2093
    9004    2093

    7103    2103
    7107    2107
    7117    2117


    % ========================================================
    % GRASSLAND / STEPPE
    % ========================================================

    7127    2126

    7135    2135
    7153    2153

    7145    2135
    7146    2135


    % ========================================================
    % RIPARIAN / WETLAND
    % ========================================================

    9011    2159
    9019    2159
    9021    2160
    9022    2180

    9327    2159
    9519    2159
    9827    2159


    % ========================================================
    % INTRODUCED / RUDERAL GRASSLAND
    % ========================================================

    9307    2181
    9308    2181
    9309    2181
    9828    2181


    % ========================================================
    % RUDERAL SHRUBLAND
    % ========================================================

    9328    2210
    9336    2210
    7078    2210


    % ========================================================
    % SAGEBRUSH / OTHER SHRUBLAND
    % ========================================================

    7125    2220
    7126    2220


    % ========================================================
    % ADDITIONAL CLASSES FOUND IN THE ACTUAL RASTER
    % ========================================================
    %
    % These were reported by the raster but were missing from
    % the previous crosswalk.
    %
    % They are assigned to the closest existing Garlick
    % habitat category used by the current model.
    %
    % ========================================================

    % --------------------------------------------------------
    % Forest / woodland
    % --------------------------------------------------------

    7056    2016
    7057    2016

    % --------------------------------------------------------
    % Shrubland
    % --------------------------------------------------------

    7070    2081

    % --------------------------------------------------------
    % Grassland / steppe
    % --------------------------------------------------------

    7143    2135
    7144    2135

    % --------------------------------------------------------
    % Developed / urban vegetation
    % --------------------------------------------------------

    7900    21
    7901    21
    7902    21

    % --------------------------------------------------------
    % Developed / agricultural / disturbed
    % --------------------------------------------------------

    7920    21
    7921    21
    7922    21
    7923    21
    7924    21

    % --------------------------------------------------------
    % Agriculture
    % --------------------------------------------------------

    7964    82

    % --------------------------------------------------------
    % Barren / disturbed
    % --------------------------------------------------------

    9017    31

    % --------------------------------------------------------
    % Ruderal / introduced vegetation
    % --------------------------------------------------------

    9329    2181
    9829    2181

];


%% ============================================================
% 5. CHECK EVERY RAW LANDFIRE CLASS
% ============================================================

fprintf('\nChecking LANDFIRE EVT classes...\n');


rawClasses = ...
    unique(EVT(EVT > 0));


mappedClasses = ...
    crosswalk(:,1);


unmatchedRaw = ...
    setdiff(rawClasses,mappedClasses);


if ~isempty(unmatchedRaw)

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' UNMATCHED LANDFIRE EVT CLASSES\n');
    fprintf('============================================\n');

    fprintf('\nThe following raw LANDFIRE EVT classes are not\n');
    fprintf('present in the crosswalk:\n\n');

    fprintf('%s\n', ...
        mat2str(unmatchedRaw'));

    error([ ...
        '\nUnmatched LANDFIRE EVT classes detected.\n' ...
        'Add these classes to the crosswalk before continuing.']);

end


fprintf('All raw LANDFIRE EVT classes are mapped.\n');


%% ============================================================
% 6. CONVERT 30-m LANDFIRE -> 30-m GARLICK EVT
% ============================================================

fprintf('\nApplying LANDFIRE -> Garlick crosswalk...\n');


PaperEVT30 = ...
    zeros(size(EVT));


for k = 1:size(crosswalk,1)

    rawCode = ...
        crosswalk(k,1);

    paperCode = ...
        crosswalk(k,2);


    mask = ...
        EVT == rawCode;


    PaperEVT30(mask) = ...
        paperCode;

end


fprintf('Crosswalk applied successfully.\n');


%% ============================================================
% 7. AGGREGATE CONVERTED LANDSCAPE TO 600 m
% ============================================================
%
% IMPORTANT:
%
% Crosswalk conversion occurs BEFORE aggregation.
%
% Each 600-m cell contains:
%
%       20 x 20 = 400
%
% 30-m pixels.
%
% The dominant Garlick EVT class is assigned to the
% 600-m model cell.
%
% ============================================================

[nRows,nCols] = ...
    size(PaperEVT30);


nModelRows = ...
    floor(nRows / aggregationFactor);


nModelCols = ...
    floor(nCols / aggregationFactor);


usedRows = ...
    nModelRows * aggregationFactor;


usedCols = ...
    nModelCols * aggregationFactor;


PaperEVT30 = ...
    PaperEVT30( ...
        1:usedRows, ...
        1:usedCols);


Landscape = ...
    zeros(nModelRows,nModelCols);


fprintf('\nAggregating to 600-m model grid...\n');


for i = 1:nModelRows

    rowStart = ...
        (i-1)*aggregationFactor + 1;

    rowEnd = ...
        i*aggregationFactor;


    for j = 1:nModelCols

        colStart = ...
            (j-1)*aggregationFactor + 1;

        colEnd = ...
            j*aggregationFactor;


        block = ...
            PaperEVT30( ...
                rowStart:rowEnd, ...
                colStart:colEnd);


        block = ...
            block(block > 0);


        if isempty(block)

            Landscape(i,j) = 0;

        else

            Landscape(i,j) = mode(block);

        end

    end

end


fprintf('600-m Garlick landscape created.\n');


fprintf('  Model grid size: %d x %d\n', ...
    nModelRows,nModelCols);

fprintf('  Physical extent: %.2f km x %.2f km\n', ...
    nModelCols * modelResolution, ...
    nModelRows * modelResolution);


%% ============================================================
% 8. GARLICK EVT MOTILITY TABLE
% ============================================================

paperEVT = [ ...
    11
    21
    22
    23
    31
    81
    82
    2001
    2006
    2011
    2016
    2051
    2052
    2054
    2055
    2057
    2061
    2062
    2064
    2066
    2080
    2081
    2086
    2093
    2103
    2107
    2117
    2126
    2135
    2153
    2159
    2160
    2180
    2181
    2210
    2214
    2217
    2220
];


%% ============================================================
% 9. FEMALE WINTER MOTILITY
% ============================================================

paperFemaleWinter = [ ...
    0.22
    0.44
    0.38
    0.41
    0.37
    0.31
    0.62
    0.51
    1.56
    0.23
    0.19
    0.20
    0.64
    0.23
    0.97
    0.18
    0.40
    0.28
    0.23
    0.31
    0.15
    0.21
    0.36
    0.64
    0.31
    0.31
    0.33
    0.12
    0.32
    0.20
    0.18
    0.12
    0.22
    0.17
    0.20
    0.22
    0.20
    0.14
];


%% ============================================================
% 10. MALE WINTER MOTILITY
% ============================================================

paperMaleWinter = [ ...
    0.10
    0.57
    0.54
    0.42
    0.44
    0.51
    1.43
    0.41
    0.82
    0.29
    0.27
    0.23
    0.53
    0.23
    1.73
    12.96
    0.66
    0.44
    0.29
    0.73
    0.20
    0.32
    0.47
    0.72
    0.13
    0.81
    0.28
    0.27
    0.28
    0.33
    0.29
    0.52
    0.34
    0.27
    0.26
    0.28
    0.21
    0.22
];


%% ============================================================
% 11. CHECK MOTILITY TABLE
% ============================================================

if length(paperEVT) ~= length(paperFemaleWinter) || ...
   length(paperEVT) ~= length(paperMaleWinter)

    error( ...
        'paperEVT and motility arrays must have equal length.');

end


%% ============================================================
% 12. CREATE WINTER MOTILITY FIELDS
% ============================================================

MuFemaleWinter = ...
    nan(nModelRows,nModelCols);


MuMaleWinter = ...
    nan(nModelRows,nModelCols);


fprintf('\nAssigning winter motility coefficients...\n');


for k = 1:length(paperEVT)

    code = ...
        paperEVT(k);


    mask = ...
        Landscape == code;


    MuFemaleWinter(mask) = ...
        paperFemaleWinter(k);


    MuMaleWinter(mask) = ...
        paperMaleWinter(k);

end


%% ============================================================
% 13. CHECK FOR PAPER CLASSES WITHOUT MOTILITY
% ============================================================

unknownMotility = ...
    (Landscape > 0) & ...
    (isnan(MuFemaleWinter) | isnan(MuMaleWinter));


if any(unknownMotility(:))

    missingCodes = ...
        unique(Landscape(unknownMotility));


    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' GARLICK EVT CLASSES WITHOUT MOTILITY\n');
    fprintf('============================================\n');

    fprintf('\n');

    disp(missingCodes');


    error([ ...
        'The crosswalk produced Garlick EVT classes ' ...
        'that are not in the motility table.']);

end


%% ============================================================
% 14. CHECK MOTILITY VALUES
% ============================================================

if any(~isfinite(MuFemaleWinter(:))) || ...
   any(MuFemaleWinter(:) <= 0)

    error( ...
        'Invalid female winter motility values detected.');

end


if any(~isfinite(MuMaleWinter(:))) || ...
   any(MuMaleWinter(:) <= 0)

    error( ...
        'Invalid male winter motility values detected.');

end


%% ============================================================
% 15. SEASONAL MOTILITY
% ============================================================
%
% At this point the available table contains the winter
% motility values.
%
% Therefore summer and breeding use the same spatial
% motility structure until their separate Table 2 values
% are inserted.
%
% The seasonal fields are still created so that the
% simulation engine can switch seasons without changing
% its interface.
%
% ============================================================

MuFemaleSummer = ...
    MuFemaleWinter;

MuFemaleBreeding = ...
    MuFemaleWinter;


MuMaleSummer = ...
    MuMaleWinter;

MuMaleBreeding = ...
    MuMaleWinter;


%% ============================================================
% 16. LOCAL 3-km HOMOGENIZATION
% ============================================================
%
% 3 km / 0.6 km = 5
%
% Therefore each local homogenization region contains
% at most 5 x 5 model cells.
%
% The harmonic mean is calculated independently for each
% local block.
%
% ============================================================

fprintf('\nPerforming local 3-km homogenization...\n');


mubarFemaleSummer = ...
    localHarmonicMean( ...
        MuFemaleSummer, ...
        homogenizationFactor);


mubarFemaleBreeding = ...
    localHarmonicMean( ...
        MuFemaleBreeding, ...
        homogenizationFactor);


mubarFemaleWinter = ...
    localHarmonicMean( ...
        MuFemaleWinter, ...
        homogenizationFactor);


mubarMaleSummer = ...
    localHarmonicMean( ...
        MuMaleSummer, ...
        homogenizationFactor);


mubarMaleBreeding = ...
    localHarmonicMean( ...
        MuMaleBreeding, ...
        homogenizationFactor);


mubarMaleWinter = ...
    localHarmonicMean( ...
        MuMaleWinter, ...
        homogenizationFactor);


fprintf('Local 3-km homogenization complete.\n');


%% ============================================================
% 17. STORE PHYSICAL MOTILITY IN PARSETS
% ============================================================

ParSets = struct();


ParSets.MuFemaleSummer = ...
    MuFemaleSummer;

ParSets.MuFemaleBreeding = ...
    MuFemaleBreeding;

ParSets.MuFemaleWinter = ...
    MuFemaleWinter;


ParSets.MuMaleSummer = ...
    MuMaleSummer;

ParSets.MuMaleBreeding = ...
    MuMaleBreeding;

ParSets.MuMaleWinter = ...
    MuMaleWinter;


%% ============================================================
% 18. STORE HOMOGENIZED MOTILITY
% ============================================================

ParSets.MuFemaleHomSummer = ...
    mubarFemaleSummer;

ParSets.MuFemaleHomBreeding = ...
    mubarFemaleBreeding;

ParSets.MuFemaleHomWinter = ...
    mubarFemaleWinter;


ParSets.MuMaleHomSummer = ...
    mubarMaleSummer;

ParSets.MuMaleHomBreeding = ...
    mubarMaleBreeding;

ParSets.MuMaleHomWinter = ...
    mubarMaleWinter;


%% ============================================================
% 19. STORE MODEL INFORMATION
% ============================================================

ParSets.modelResolution = ...
    modelResolution;

ParSets.homogenizationResolution = ...
    homogenizationResolution;

ParSets.homogenizationFactor = ...
    homogenizationFactor;

ParSets.Landscape = ...
    Landscape;


%% ============================================================
% 20. CREATE MU STRUCTURE
% ============================================================

Mu = struct();


Mu.FemaleSummer = ...
    MuFemaleSummer;

Mu.FemaleBreeding = ...
    MuFemaleBreeding;

Mu.FemaleWinter = ...
    MuFemaleWinter;


Mu.MaleSummer = ...
    MuMaleSummer;

Mu.MaleBreeding = ...
    MuMaleBreeding;

Mu.MaleWinter = ...
    MuMaleWinter;


Mu.FemaleHomSummer = ...
    mubarFemaleSummer;

Mu.FemaleHomBreeding = ...
    mubarFemaleBreeding;

Mu.FemaleHomWinter = ...
    mubarFemaleWinter;


Mu.MaleHomSummer = ...
    mubarMaleSummer;

Mu.MaleHomBreeding = ...
    mubarMaleBreeding;

Mu.MaleHomWinter = ...
    mubarMaleWinter;


%% ============================================================
% 21. CREATE MUBAR STRUCTURE
% ============================================================

mubar = struct();


mubar.FemaleSummer = ...
    mubarFemaleSummer;

mubar.FemaleBreeding = ...
    mubarFemaleBreeding;

mubar.FemaleWinter = ...
    mubarFemaleWinter;


mubar.MaleSummer = ...
    mubarMaleSummer;

mubar.MaleBreeding = ...
    mubarMaleBreeding;

mubar.MaleWinter = ...
    mubarMaleWinter;


%% ============================================================
% 22. INITIAL SUMMER SUSCEPTIBLE POPULATIONS
% ============================================================
%
% Main begins in summer.
%
%       SF = CF / MuFemaleSummer
%
%       SM = CM / MuMaleSummer
%
% CF and CM are selected so that:
%
%       integral(SF) = NF
%
%       integral(SM) = NM
%
% ============================================================

invMuF = ...
    1 ./ MuFemaleSummer;


invMuM = ...
    1 ./ MuMaleSummer;


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


SS0A = ...
    SF + SM;


%% ============================================================
% 23. CHECK POPULATIONS
% ============================================================

NF_check = ...
    sum(SF(:)) * cellArea;


NM_check = ...
    sum(SM(:)) * cellArea;


fprintf('\n');
fprintf('============================================\n');
fprintf(' INITIAL GARLICK POPULATION\n');
fprintf('============================================\n');


fprintf('\nFemale population:\n');

fprintf('  Requested:  %.4f\n', ...
    NF);

fprintf('  Calculated: %.4f\n', ...
    NF_check);


fprintf('\nMale population:\n');

fprintf('  Requested:  %.4f\n', ...
    NM);

fprintf('  Calculated: %.4f\n', ...
    NM_check);


fprintf('\nTotal population:\n');

fprintf('  Requested:  %.4f\n', ...
    NF + NM);

fprintf('  Calculated: %.4f\n', ...
    NF_check + NM_check);


%% ============================================================
% 24. MOTILITY SUMMARY
% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf(' MOTILITY SUMMARY\n');
fprintf('============================================\n');


fprintf('\nFemale:\n');

fprintf('  Summer:   %.6g - %.6g\n', ...
    min(MuFemaleSummer(:)), ...
    max(MuFemaleSummer(:)));

fprintf('  Breeding: %.6g - %.6g\n', ...
    min(MuFemaleBreeding(:)), ...
    max(MuFemaleBreeding(:)));

fprintf('  Winter:   %.6g - %.6g\n', ...
    min(MuFemaleWinter(:)), ...
    max(MuFemaleWinter(:)));


fprintf('\nMale:\n');

fprintf('  Summer:   %.6g - %.6g\n', ...
    min(MuMaleSummer(:)), ...
    max(MuMaleSummer(:)));

fprintf('  Breeding: %.6g - %.6g\n', ...
    min(MuMaleBreeding(:)), ...
    max(MuMaleBreeding(:)));

fprintf('  Winter:   %.6g - %.6g\n', ...
    min(MuMaleWinter(:)), ...
    max(MuMaleWinter(:)));


fprintf('\n');
fprintf('Spatial background complete.\n');


%% ============================================================
% CLOSE MAIN FUNCTION
% ============================================================

end


%% ============================================================
% LOCAL FUNCTION
% LOCAL 3-km HARMONIC MEAN
% ============================================================

function mubar = ...
    localHarmonicMean(mu,blockSize)


[nRows,nCols] = ...
    size(mu);


mubar = ...
    zeros(nRows,nCols);


for rStart = 1:blockSize:nRows

    rEnd = ...
        min(rStart + blockSize - 1,nRows);


    for cStart = 1:blockSize:nCols

        cEnd = ...
            min(cStart + blockSize - 1,nCols);


        block = ...
            mu(rStart:rEnd,cStart:cEnd);


        valid = ...
            isfinite(block) & ...
            block > 0;


        values = ...
            block(valid);


        if isempty(values)

            mubar( ...
                rStart:rEnd, ...
                cStart:cEnd) = NaN;

        else

            harmonicValue = ...
                numel(values) / ...
                sum(1 ./ values);


            mubar( ...
                rStart:rEnd, ...
                cStart:cEnd) = ...
                harmonicValue;

        end

    end

end

end