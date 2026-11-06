function results = CWDSimulationGarlick( ...
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
    CM)

% ============================================================
% HOMOGENIZED GARLICK et al. (2014) CWD MODEL
% ============================================================
%
% This implementation uses the homogenized formulation of
% Garlick et al. (2014), but keeps the simulation in PHYSICAL
% TIME (days).
%
% IMPORTANT:
%
% The equations printed as Eqs. 39-41 in Garlick et al. contain
% 1/omega factors because they are written in the rescaled
% homogenized formulation.
%
% Those 1/omega factors should NOT simply be inserted into a
% physical 5-day timestep. Doing so makes the numerical system
% extremely stiff and produces Inf/NaN values.
%
% This file therefore uses the corresponding physical-time
% formulation.
%
%
% STATE VARIABLES
%
%   cIF = homogenized infected-female variable
%   cIM = homogenized infected-male variable
%   H   = environmental hazard
%
%
% PHYSICAL INFECTED DENSITIES
%
%   IF = cIF ./ MuFemale
%   IM = cIM ./ MuMale
%
%
% SUSCEPTIBLE STEADY-STATE APPROXIMATION
%
%   SF = CF ./ MuFemale
%   SM = CM ./ MuMale
%
%
% The susceptible populations are held at the Garlick
% steady-state values supplied by the main program.
%
%
% NUMERICAL METHOD
%
%   Reaction: exact matrix exponential
%   Diffusion: ADI
%   Splitting: reaction half-step
%               diffusion full-step
%               reaction half-step
%
% ============================================================


%% ============================================================
% 1. PARAMETERS
% ============================================================

dt = params.dt;
T  = params.T;

dx = params.dx;
dy = params.dy;

cellArea = params.cellArea;


betaFF = params.betaFF;
betaMF = params.betaMF;

betaFM = params.betaFM;
betaMM = params.betaMM;


gammaF = params.gammaF;
gammaM = params.gammaM;


alphaF = params.alphaF;
alphaM = params.alphaM;

alphaDF = params.alphaDF;
alphaDM = params.alphaDM;


omega = params.omega;

delta = params.delta;

phi = params.phi;


%% ============================================================
% 2. CURRENT MOTILITY
% ============================================================

MuFemale = params.MuFemale;
MuMale   = params.MuMale;


% Make sure the motility fields are valid.

if any(~isfinite(MuFemale(:))) || ...
   any(MuFemale(:) <= 0)

    error( ...
        'Female motility contains zero, negative, or nonfinite values.');

end


if any(~isfinite(MuMale(:))) || ...
   any(MuMale(:) <= 0)

    error( ...
        'Male motility contains zero, negative, or nonfinite values.');

end


%% ============================================================
% 3. HOMOGENIZED MOTILITIES
% ============================================================

mubarF = getMubar(mubar,'Female');

mubarM = getMubar(mubar,'Male');


% If buildBackground did not provide scalar values,
% calculate harmonic means from the current motility fields.

if isempty(mubarF)

    mubarF = ...
        1 / mean(1 ./ MuFemale(:));

end


if isempty(mubarM)

    mubarM = ...
        1 / mean(1 ./ MuMale(:));

end


% If the supplied value is a spatial field, reduce it to
% a scalar homogenized coefficient.

if ~isscalar(mubarF)

    mubarF = ...
        mean(mubarF(:),'omitnan');

end


if ~isscalar(mubarM)

    mubarM = ...
        mean(mubarM(:),'omitnan');

end


if ~isfinite(mubarF) || mubarF <= 0

    error('Invalid female homogenized motility.');

end


if ~isfinite(mubarM) || mubarM <= 0

    error('Invalid male homogenized motility.');

end


fprintf('\n');
fprintf('============================================\n');
fprintf(' HOMOGENIZED MOTILITIES\n');
fprintf('============================================\n');

fprintf('\nFemale mubar: %.8g km^2/day\n',mubarF);
fprintf('Male mubar:   %.8g km^2/day\n',mubarM);


%% ============================================================
% 4. STEADY-STATE SUSCEPTIBLES
% ============================================================

SF = SF0;
SM = SM0;


checkSF = ...
    sum(SF(:)) * cellArea;

checkSM = ...
    sum(SM(:)) * cellArea;


fprintf('\nSteady-state female population: %.4f\n',checkSF);
fprintf('Steady-state male population:   %.4f\n',checkSM);


%% ============================================================
% 5. INITIAL INFECTED VARIABLES
% ============================================================

cIF = MuFemale .* IF0;

cIM = MuMale .* IM0;

H = H0;


% Remove tiny numerical negative values.

cIF = max(real(cIF),0);
cIM = max(real(cIM),0);
H   = max(real(H),0);


%% ============================================================
% 6. INITIAL POPULATION
% ============================================================

initialIF = ...
    sum(IF0(:)) * cellArea;

initialIM = ...
    sum(IM0(:)) * cellArea;


%% ============================================================
% 7. SIMULATION SETTINGS
% ============================================================

N = round(T/dt);

outputDt = params.outputDt;

saveEvery = ...
    max(1,round(outputDt/dt));

numSaved = ...
    floor(N/saveEvery) + 1;


%% ============================================================
% 8. GRID SIZE
% ============================================================

[L1,L2] = size(SF);


if ~isequal(size(SF),size(SM), ...
            size(cIF),size(cIM),size(H))

    error( ...
        'SF, SM, cIF, cIM, and H must have identical dimensions.');

end


%% ============================================================
% 9. PREALLOCATE OUTPUTS
% ============================================================

time = zeros(numSaved,1);


IFsave = ...
    zeros(L1,L2,numSaved);

IMsave = ...
    zeros(L1,L2,numSaved);

Hsave = ...
    zeros(L1,L2,numSaved);


SFsave = ...
    zeros(L1,L2,numSaved);

SMsave = ...
    zeros(L1,L2,numSaved);


totalIF = zeros(numSaved,1);
totalIM = zeros(numSaved,1);

totalSF = zeros(numSaved,1);
totalSM = zeros(numSaved,1);

totalS = zeros(numSaved,1);
totalI = zeros(numSaved,1);

totalDeer = zeros(numSaved,1);

totalH = zeros(numSaved,1);

prevalence = zeros(numSaved,1);


%% ============================================================
% 10. GARLICK SUSCEPTIBLE DENSITIES
% ============================================================

S_hatF = ...
    CF ./ MuFemale;

S_hatM = ...
    CM ./ MuMale;


if any(~isfinite(S_hatF(:))) || ...
   any(~isfinite(S_hatM(:)))

    error( ...
        ['Nonfinite susceptible steady-state density detected. ' ...
         'Check MuFemale, MuMale, CF, and CM.']);

end


%% ============================================================
% 11. PHYSICAL-TIME HOMOGENIZED REACTION TERMS
% ============================================================
%
% IMPORTANT:
%
% We are NOT using the raw 1/omega-scaled coefficients here.
%
% Instead we retain the homogenized structure while keeping
% time in physical days.
%
%
% Female infected:
%
% dcIF/dt =
%
%   mubarF^2 * Laplacian(cIF)
%
%   + A_F*cIF
%   + B_FM*cIM
%   + E_F*H
%
%
% Male infected:
%
% dcIM/dt =
%
%   mubarM^2 * Laplacian(cIM)
%
%   + B_MF*cIF
%   + A_M*cIM
%   + E_M*H
%
%
% Hazard:
%
% dH/dt =
%
%   Q_F*cIF
%   + Q_M*cIM
%   - delta*H
%
%
% ============================================================


%% ------------------------------------------------------------
% Female infected self-transmission
% ------------------------------------------------------------

A_F = ...
    mubarF .* ...
    mean( ...
        (betaFF .* S_hatF - omega - phi) ...
        ./ MuFemale, ...
        'all');


%% ------------------------------------------------------------
% Male infected -> female infected
% ------------------------------------------------------------

B_FM = ...
    mubarF .* ...
    mean( ...
        betaMF .* S_hatF ...
        ./ MuMale, ...
        'all');


%% ------------------------------------------------------------
% Environmental -> female infection
% ------------------------------------------------------------

E_F = ...
    mubarF .* ...
    mean( ...
        gammaF .* alphaDM .* S_hatF, ...
        'all');


%% ------------------------------------------------------------
% Male infected self-transmission
% ------------------------------------------------------------

A_M = ...
    mubarM .* ...
    mean( ...
        (betaMM .* S_hatM - omega - phi) ...
        ./ MuMale, ...
        'all');


%% ------------------------------------------------------------
% Female infected -> male infected
% ------------------------------------------------------------

B_MF = ...
    mubarM .* ...
    mean( ...
        betaFM .* S_hatM ...
        ./ MuFemale, ...
        'all');


%% ------------------------------------------------------------
% Environmental -> male infection
% ------------------------------------------------------------

E_M = ...
    mubarM .* ...
    mean( ...
        gammaM .* alphaDM .* S_hatM, ...
        'all');


%% ------------------------------------------------------------
% Environmental contamination from females
% ------------------------------------------------------------

meanInvMuF = ...
    mean(1 ./ MuFemale(:));


meanInvMuM = ...
    mean(1 ./ MuMale(:));


Q_F = ...
    meanInvMuF .* ...
    (alphaF + alphaDF .* omega);


%% ------------------------------------------------------------
% Environmental contamination from males
% ------------------------------------------------------------

Q_M = ...
    meanInvMuM .* ...
    (alphaM + alphaDM .* omega);


%% ------------------------------------------------------------
% Environmental decay
% ------------------------------------------------------------

Q_H = -delta;


%% ============================================================
% 12. DISPLAY COEFFICIENTS
% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf(' PHYSICAL-TIME REACTION COEFFICIENTS\n');
fprintf('============================================\n');

fprintf('\nA_F  = %.8g\n',A_F);
fprintf('B_FM = %.8g\n',B_FM);
fprintf('E_F  = %.8g\n',E_F);

fprintf('\nA_M  = %.8g\n',A_M);
fprintf('B_MF = %.8g\n',B_MF);
fprintf('E_M  = %.8g\n',E_M);

fprintf('\nQ_F  = %.8g\n',Q_F);
fprintf('Q_M  = %.8g\n',Q_M);
fprintf('Q_H  = %.8g\n',Q_H);


%% ============================================================
% 13. REACTION MATRIX
% ============================================================
%
% The reaction system is linear:
%
%       d/dt [cIF] = A [cIF]
%          [cIM]       [cIM]
%          [ H ]       [ H ]
%
% Because the coefficients are spatially homogenized,
% the same 3x3 matrix is applied at every model cell.
%
% ============================================================

A = [ ...
    A_F,  B_FM, E_F;
    B_MF, A_M,  E_M;
    Q_F,  Q_M,  Q_H];


%% ============================================================
% 14. CHECK REACTION MATRIX
% ============================================================

if any(~isfinite(A(:)))

    error( ...
        'Reaction matrix contains Inf or NaN.');

end


fprintf('\nReaction matrix:\n');

disp(A);


eigA = eig(A);


fprintf('Reaction eigenvalues:\n');

for k = 1:length(eigA)

    fprintf('  %.8g\n',eigA(k));

end


%% ============================================================
% 15. DIFFUSION COEFFICIENTS
% ============================================================
%
% Physical-time homogenized diffusion:
%
%       D_F = mubarF^2
%       D_M = mubarM^2
%
% NOT:
%
%       mubar^2 / omega
%
% The latter belongs to the rescaled equation and caused the
% numerical explosion when combined with physical days.
%
% ============================================================

D_F = ...
    mubarF^2;

D_M = ...
    mubarM^2;


fprintf('\n');
fprintf('============================================\n');
fprintf(' DIFFUSION\n');
fprintf('============================================\n');

fprintf('\nFemale D = %.8g\n',D_F);
fprintf('Male D   = %.8g\n',D_M);


%% ============================================================
% 16. ADI OPERATORS
% ============================================================

AxF = ...
    buildADIoperator( ...
        L2,dx,D_F,dt);

BxF = ...
    buildADIexplicit( ...
        L2,dx,D_F,dt);

AyF = ...
    buildADIoperator( ...
        L1,dy,D_F,dt);

ByF = ...
    buildADIexplicit( ...
        L1,dy,D_F,dt);


AxM = ...
    buildADIoperator( ...
        L2,dx,D_M,dt);

BxM = ...
    buildADIexplicit( ...
        L2,dx,D_M,dt);

AyM = ...
    buildADIoperator( ...
        L1,dy,D_M,dt);

ByM = ...
    buildADIexplicit( ...
        L1,dy,D_M,dt);


%% ============================================================
% 17. EXACT REACTION OPERATOR
% ============================================================
%
% This is the major numerical fix.
%
% Instead of RK4 on the very stiff rescaled equations, use
% the exact matrix exponential of the physical-time reaction
% matrix.
%
% ============================================================

reactionHalf = ...
    expm(A * (dt/2));


if any(~isfinite(reactionHalf(:)))

    error( ...
        'Reaction matrix exponential is nonfinite.');

end


%% ============================================================
% 18. SAVE INITIAL STATE
% ============================================================

saveIndex = 1;


time(saveIndex) = 0;


IFsave(:,:,saveIndex) = ...
    cIF ./ MuFemale;

IMsave(:,:,saveIndex) = ...
    cIM ./ MuMale;

Hsave(:,:,saveIndex) = H;


SFsave(:,:,saveIndex) = SF;
SMsave(:,:,saveIndex) = SM;


[ ...
    totalSF(saveIndex), ...
    totalSM(saveIndex), ...
    totalIF(saveIndex), ...
    totalIM(saveIndex), ...
    totalS(saveIndex), ...
    totalI(saveIndex), ...
    totalDeer(saveIndex), ...
    totalH(saveIndex), ...
    prevalence(saveIndex)] = ...
    calculateGarlickTotals( ...
        SF, ...
        SM, ...
        cIF, ...
        cIM, ...
        H, ...
        MuFemale, ...
        MuMale, ...
        cellArea);


%% ============================================================
% 19. SIMULATION START
% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf(' STARTING GARLICK SIMULATION\n');
fprintf('============================================\n');

fprintf('\nSimulation length: %.2f days\n',T);
fprintf('Time step: %.2f days\n',dt);

fprintf('\nInitial infected females: %.4f\n',initialIF);
fprintf('Initial infected males:   %.4f\n',initialIM);

fprintf('\n');


%% ============================================================
% 20. SIMULATION LOOP
% ============================================================

for n = 1:N


    %% --------------------------------------------------------
    % REACTION HALF STEP
    % --------------------------------------------------------

    [cIF,cIM,H] = ...
        applyReaction( ...
            cIF, ...
            cIM, ...
            H, ...
            reactionHalf);


    %% --------------------------------------------------------
    % DIFFUSION
    % --------------------------------------------------------

    cIF = ...
        diffuseADI( ...
            cIF, ...
            AxF, ...
            BxF, ...
            AyF, ...
            ByF);


    cIM = ...
        diffuseADI( ...
            cIM, ...
            AxM, ...
            BxM, ...
            AyM, ...
            ByM);


    %% --------------------------------------------------------
    % REACTION HALF STEP
    % --------------------------------------------------------

    [cIF,cIM,H] = ...
        applyReaction( ...
            cIF, ...
            cIM, ...
            H, ...
            reactionHalf);


    %% --------------------------------------------------------
    % POSITIVITY
    % --------------------------------------------------------

    cIF = max(real(cIF),0);
    cIM = max(real(cIM),0);
    H   = max(real(H),0);


    %% --------------------------------------------------------
    % NUMERICAL CHECK
    % --------------------------------------------------------

    if any(~isfinite(cIF(:)))

        error( ...
            ['Female infection became Inf/NaN at ' ...
             't = %.2f days.'], ...
             n*dt);

    end


    if any(~isfinite(cIM(:)))

        error( ...
            ['Male infection became Inf/NaN at ' ...
             't = %.2f days.'], ...
             n*dt);

    end


    if any(~isfinite(H(:)))

        error( ...
            ['Environmental hazard became Inf/NaN at ' ...
             't = %.2f days.'], ...
             n*dt);

    end


    %% --------------------------------------------------------
    % SAVE OUTPUT
    % --------------------------------------------------------

    if mod(n,saveEvery) == 0 || n == N

        saveIndex = saveIndex + 1;


        time(saveIndex) = ...
            n*dt;


        IFsave(:,:,saveIndex) = ...
            cIF ./ MuFemale;

        IMsave(:,:,saveIndex) = ...
            cIM ./ MuMale;

        Hsave(:,:,saveIndex) = ...
            H;


        SFsave(:,:,saveIndex) = ...
            SF;

        SMsave(:,:,saveIndex) = ...
            SM;


        [ ...
            totalSF(saveIndex), ...
            totalSM(saveIndex), ...
            totalIF(saveIndex), ...
            totalIM(saveIndex), ...
            totalS(saveIndex), ...
            totalI(saveIndex), ...
            totalDeer(saveIndex), ...
            totalH(saveIndex), ...
            prevalence(saveIndex)] = ...
            calculateGarlickTotals( ...
                SF, ...
                SM, ...
                cIF, ...
                cIM, ...
                H, ...
                MuFemale, ...
                MuMale, ...
                cellArea);

    end


    %% --------------------------------------------------------
    % PROGRESS
    % --------------------------------------------------------

    if mod(n,max(1,round(N/10))) == 0

        fprintf( ...
            ' %5.1f%% complete\n', ...
            100*n/N);

    end

end


%% ============================================================
% 21. TRIM OUTPUTS
% ============================================================

time = ...
    time(1:saveIndex);


SFsave = ...
    SFsave(:,:,1:saveIndex);

SMsave = ...
    SMsave(:,:,1:saveIndex);


IFsave = ...
    IFsave(:,:,1:saveIndex);

IMsave = ...
    IMsave(:,:,1:saveIndex);

Hsave = ...
    Hsave(:,:,1:saveIndex);


totalSF = ...
    totalSF(1:saveIndex);

totalSM = ...
    totalSM(1:saveIndex);

totalIF = ...
    totalIF(1:saveIndex);

totalIM = ...
    totalIM(1:saveIndex);

totalS = ...
    totalS(1:saveIndex);

totalI = ...
    totalI(1:saveIndex);

totalDeer = ...
    totalDeer(1:saveIndex);

totalH = ...
    totalH(1:saveIndex);

prevalence = ...
    prevalence(1:saveIndex);


%% ============================================================
% 22. RESULTS STRUCTURE
% ============================================================

results = struct();


results.time = time;


results.SF = SFsave;
results.SM = SMsave;


results.IF = IFsave;
results.IM = IMsave;


results.H = Hsave;


results.cIF = ...
    IFsave .* MuFemale;

results.cIM = ...
    IMsave .* MuMale;


results.totalSF = totalSF;
results.totalSM = totalSM;


results.totalIF = totalIF;
results.totalIM = totalIM;


results.totalS = totalS;
results.totalI = totalI;


results.totalDeer = totalDeer;

results.totalH = totalH;

results.prevalence = prevalence;


results.MuFemale = MuFemale;
results.MuMale = MuMale;


results.mubarFemale = mubarF;
results.mubarMale = mubarM;


results.Mu = Mu;
results.mubar = mubar;


results.Landscape = Landscape;


results.cellArea = cellArea;

results.dx = dx;
results.dy = dy;

results.outputDt = outputDt;


results.CF = CF;
results.CM = CM;


%% ============================================================
% 23. FINAL SUMMARY
% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf(' SIMULATION COMPLETE\n');
fprintf('============================================\n');

fprintf('\nFinal populations:\n');

fprintf(' Susceptible females: %.4f\n', ...
    totalSF(end));

fprintf(' Susceptible males:   %.4f\n', ...
    totalSM(end));

fprintf(' Infected females:    %.4f\n', ...
    totalIF(end));

fprintf(' Infected males:      %.4f\n', ...
    totalIM(end));

fprintf('\n Total infected:      %.4f\n', ...
    totalI(end));

fprintf(' Total deer:          %.4f\n', ...
    totalDeer(end));

fprintf(' Final prevalence:    %.6f\n', ...
    prevalence(end));


end


%% ============================================================
% APPLY EXACT REACTION STEP
% ============================================================

function [cIF,cIM,H] = ...
    applyReaction(cIF,cIM,H,E)

% ------------------------------------------------------------
% Apply the same 3x3 reaction matrix independently at every
% spatial location.
%
% E is the matrix exponential:
%
%       E = expm(A*dt/2)
%
% ------------------------------------------------------------

Fnew = ...
      E(1,1).*cIF ...
    + E(1,2).*cIM ...
    + E(1,3).*H;


Mnew = ...
      E(2,1).*cIF ...
    + E(2,2).*cIM ...
    + E(2,3).*H;


Hnew = ...
      E(3,1).*cIF ...
    + E(3,2).*cIM ...
    + E(3,3).*H;


cIF = Fnew;
cIM = Mnew;
H   = Hnew;


end


%% ============================================================
% GET MUBAR
% ============================================================

function value = getMubar(mubar,sex)

value = [];


if isempty(mubar)

    return;

end


if isnumeric(mubar)

    if isscalar(mubar)

        value = mubar;

    end

    return;

end


if ~isstruct(mubar)

    return;

end


if strcmpi(sex,'Female')

    names = { ...
        'Female', ...
        'MuFemale', ...
        'mubarFemale', ...
        'mubarF', ...
        'FemaleSummer', ...
        'MuFemaleSummer'};

else

    names = { ...
        'Male', ...
        'MuMale', ...
        'mubarMale', ...
        'mubarM', ...
        'MaleSummer', ...
        'MuMaleSummer'};

end


for i = 1:length(names)

    if isfield(mubar,names{i})

        value = ...
            mubar.(names{i});

        return;

    end

end


end


%% ============================================================
% BUILD ADI IMPLICIT OPERATOR
% ============================================================

function A = ...
    buildADIoperator(n,h,D,dt)

L = ...
    buildNeumannLaplacian(n,h);


A = ...
    speye(n) ...
    - (dt/2)*D*L;


end


%% ============================================================
% BUILD ADI EXPLICIT OPERATOR
% ============================================================

function B = ...
    buildADIexplicit(n,h,D,dt)

L = ...
    buildNeumannLaplacian(n,h);


B = ...
    speye(n) ...
    + (dt/2)*D*L;


end


%% ============================================================
% 1-D NEUMANN LAPLACIAN
% ============================================================

function L = ...
    buildNeumannLaplacian(n,h)

if n == 1

    L = sparse(1,1);

    return;

end


L = sparse(n,n);


if n > 2

    idx = 2:n-1;


    L(sub2ind([n,n],idx,idx-1)) = 1;

    L(sub2ind([n,n],idx,idx)) = -2;

    L(sub2ind([n,n],idx,idx+1)) = 1;

end


% Left boundary

L(1,1) = -2;

L(1,2) = 2;


% Right boundary

L(n,n) = -2;

L(n,n-1) = 2;


L = ...
    L / h^2;


end


%% ============================================================
% ADI DIFFUSION
% ============================================================

function Unew = ...
    diffuseADI(U,Ax,Bx,Ay,By)

% First half-step in y

rhs = ...
    By * U;


Ustar = ...
    rhs / Ax';


% Second half-step in x

rhs = ...
    Ustar * Bx';


Unew = ...
    Ay \ rhs;


Unew = ...
    real(Unew);


Unew = ...
    max(Unew,0);


end


%% ============================================================
% CALCULATE TOTALS
% ============================================================

function [ ...
    totalSF, ...
    totalSM, ...
    totalIF, ...
    totalIM, ...
    totalS, ...
    totalI, ...
    totalDeer, ...
    totalH, ...
    prevalence] = ...
    calculateGarlickTotals( ...
        SF, ...
        SM, ...
        cIF, ...
        cIM, ...
        H, ...
        MuFemale, ...
        MuMale, ...
        cellArea)


IF = ...
    cIF ./ MuFemale;


IM = ...
    cIM ./ MuMale;


totalSF = ...
    sum(SF(:))*cellArea;


totalSM = ...
    sum(SM(:))*cellArea;


totalIF = ...
    sum(IF(:))*cellArea;


totalIM = ...
    sum(IM(:))*cellArea;


totalS = ...
    totalSF + totalSM;


totalI = ...
    totalIF + totalIM;


totalDeer = ...
    totalS + totalI;


totalH = ...
    sum(H(:))*cellArea;


prevalence = ...
    totalI / (totalDeer + eps);


end