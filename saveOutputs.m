function saveOutputs(results,Landscape,tifFile)

% ============================================================
% saveOutputs
%
% Saves and plots the results from the spatial CWD SIH model.
%
% Outputs include:
%
%   1. Total deer population
%   2. Sex-specific populations
%   3. CWD prevalence
%   4. LANDFIRE landscape
%   5. Female and male motility
%   6. Final susceptible population
%   7. Final infected population
%   8. Final environmental hazard
%   9. Sex-specific final populations
%  10. Complete MATLAB results file
%
% All outputs are saved in:
%
%       outputs/
%
% ============================================================


%% ============================================================
% 1. CREATE OUTPUT DIRECTORY
% ============================================================

outputFolder = 'outputs';

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

fprintf('\nOutput directory:\n');
fprintf('  %s\n',outputFolder);


%% ============================================================
% 2. EXTRACT RESULTS
% ============================================================

time = results.time;

SF = results.SF;
SM = results.SM;

IF = results.IF;
IM = results.IM;

H = results.H;

totalSF = results.totalSF;
totalSM = results.totalSM;

totalIF = results.totalIF;
totalIM = results.totalIM;

totalS = results.totalS;
totalI = results.totalI;

totalDeer = results.totalDeer;
totalH = results.totalH;

prevalence = results.prevalence;

MuFemale = results.MuFemale;
MuMale   = results.MuMale;

Mu = results.Mu;


%% ============================================================
% 3. FINAL SPATIAL FIELDS
% ============================================================

SfinalF = SF(:,:,end);
SfinalM = SM(:,:,end);

IfinalF = IF(:,:,end);
IfinalM = IM(:,:,end);

Sfinal = SfinalF + SfinalM;
Ifinal = IfinalF + IfinalM;

Hfinal = H(:,:,end);


%% ============================================================
% 4. SAVE COMPLETE RESULTS
% ============================================================

fprintf('\nSaving complete simulation results...\n');

save( ...
    fullfile(outputFolder,'CWD_Simulation_Results.mat'), ...
    'results', ...
    '-v7.3');


%% ============================================================
% 5. SAVE SUMMARY DATA
% ============================================================

summary = struct();

summary.finalSusceptibleFemales = totalSF(end);
summary.finalSusceptibleMales   = totalSM(end);

summary.finalInfectedFemales = totalIF(end);
summary.finalInfectedMales   = totalIM(end);

summary.finalSusceptible = totalS(end);
summary.finalInfected    = totalI(end);

summary.finalTotalDeer = totalDeer(end);
summary.finalHazard    = totalH(end);

summary.finalPrevalence = prevalence(end);

summary.initialTotalDeer = totalDeer(1);
summary.initialTotalInfected = totalI(1);

summary.tifFile = tifFile;

save( ...
    fullfile(outputFolder,'Simulation_Summary.mat'), ...
    'summary');


%% ============================================================
% 6. TOTAL DEER POPULATION
% ============================================================

figure('Visible','off');

plot(time,totalDeer,'LineWidth',2);

xlabel('Time');
ylabel('Total deer');

title('Total Deer Population');

grid on;

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Total_Deer_Population.png'));

close(gcf);


%% ============================================================
% 7. SUSCEPTIBLE AND INFECTED POPULATION
% ============================================================

figure('Visible','off');

plot(time,totalS,'LineWidth',2);
hold on;

plot(time,totalI,'LineWidth',2);

xlabel('Time');
ylabel('Population');

title('Susceptible and Infected Deer');

legend( ...
    'Susceptible', ...
    'Infected', ...
    'Location','best');

grid on;

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Susceptible_Infected_Population.png'));

close(gcf);


%% ============================================================
% 8. SEX-SPECIFIC POPULATIONS
% ============================================================

figure('Visible','off');

plot(time,totalSF,'LineWidth',2);
hold on;

plot(time,totalSM,'LineWidth',2);

plot(time,totalIF,'LineWidth',2);
plot(time,totalIM,'LineWidth',2);

xlabel('Time');
ylabel('Population');

title('Sex-Specific Deer Populations');

legend( ...
    'Susceptible Females', ...
    'Susceptible Males', ...
    'Infected Females', ...
    'Infected Males', ...
    'Location','best');

grid on;

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Sex_Specific_Populations.png'));

close(gcf);


%% ============================================================
% 9. CWD PREVALENCE
% ============================================================

figure('Visible','off');

plot(time,100 .* prevalence,'LineWidth',2);

xlabel('Time');
ylabel('Prevalence (%)');

title('CWD Prevalence');

grid on;

saveas( ...
    gcf, ...
    fullfile(outputFolder,'CWD_Prevalence.png'));

close(gcf);


%% ============================================================
% 10. ENVIRONMENTAL HAZARD
% ============================================================

figure('Visible','off');

plot(time,totalH,'LineWidth',2);

xlabel('Time');
ylabel('Total environmental hazard');

title('Environmental Hazard');

grid on;

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Environmental_Hazard.png'));

close(gcf);


%% ============================================================
% 11. LANDFIRE LANDSCAPE
% ============================================================

figure('Visible','off');

imagesc(Landscape);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Aggregated LANDFIRE EVT Landscape');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'LANDFIRE_Landscape.png'));

close(gcf);


%% ============================================================
% 12. FEMALE MOTILITY
% ============================================================

figure('Visible','off');

imagesc(MuFemale);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Female Habitat-Dependent Motility');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Female_Motility.png'));

close(gcf);


%% ============================================================
% 13. MALE MOTILITY
% ============================================================

figure('Visible','off');

imagesc(MuMale);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Male Habitat-Dependent Motility');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Male_Motility.png'));

close(gcf);


%% ============================================================
% 14. FINAL SUSCEPTIBLE POPULATION
% ============================================================

figure('Visible','off');

imagesc(Sfinal);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Final Susceptible Deer');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Final_Susceptible.png'));

close(gcf);


%% ============================================================
% 15. FINAL INFECTED POPULATION
% ============================================================

figure('Visible','off');

imagesc(Ifinal);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Final Infected Deer');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Final_Infected.png'));

close(gcf);


%% ============================================================
% 16. FINAL ENVIRONMENTAL HAZARD
% ============================================================

figure('Visible','off');

imagesc(Hfinal);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Final Environmental Hazard');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Final_Hazard.png'));

close(gcf);


%% ============================================================
% 17. FINAL SUSCEPTIBLE FEMALES
% ============================================================

figure('Visible','off');

imagesc(SfinalF);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Final Susceptible Females');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Final_Susceptible_Females.png'));

close(gcf);


%% ============================================================
% 18. FINAL SUSCEPTIBLE MALES
% ============================================================

figure('Visible','off');

imagesc(SfinalM);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Final Susceptible Males');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Final_Susceptible_Males.png'));

close(gcf);


%% ============================================================
% 19. FINAL INFECTED FEMALES
% ============================================================

figure('Visible','off');

imagesc(IfinalF);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Final Infected Females');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Final_Infected_Females.png'));

close(gcf);


%% ============================================================
% 20. FINAL INFECTED MALES
% ============================================================

figure('Visible','off');

imagesc(IfinalM);

axis image;
colorbar;

xlabel('Column');
ylabel('Row');

title('Final Infected Males');

saveas( ...
    gcf, ...
    fullfile(outputFolder,'Final_Infected_Males.png'));

close(gcf);


%% ============================================================
% 21. SAVE FINAL SPATIAL DATA
% ============================================================

save( ...
    fullfile(outputFolder,'Final_Spatial_Fields.mat'), ...
    'Sfinal', ...
    'Ifinal', ...
    'Hfinal', ...
    'SfinalF', ...
    'SfinalM', ...
    'IfinalF', ...
    'IfinalM', ...
    'MuFemale', ...
    'MuMale', ...
    'Mu', ...
    'Landscape');


%% ============================================================
% 22. SAVE FINAL FIELDS AS TIFF FILES
% ============================================================

% These are standard TIFF files containing the model values.
%
% They are NOT GeoTIFF files because the current cropped
% LANDFIRE preprocessing does not preserve geographic
% referencing information.

fprintf('\nSaving final spatial TIFF files...\n');

% ------------------------------------------------------------
% Combined susceptible
% ------------------------------------------------------------

Sout = uint16(65535 .* mat2gray(Sfinal));

imwrite( ...
    Sout, ...
    fullfile(outputFolder,'Susceptible_Final.tif'));


% ------------------------------------------------------------
% Combined infected
% ------------------------------------------------------------

Iout = uint16(65535 .* mat2gray(Ifinal));

imwrite( ...
    Iout, ...
    fullfile(outputFolder,'Infected_Final.tif'));


% ------------------------------------------------------------
% Environmental hazard
% ------------------------------------------------------------

Hout = uint16(65535 .* mat2gray(Hfinal));

imwrite( ...
    Hout, ...
    fullfile(outputFolder,'Hazard_Final.tif'));


% ------------------------------------------------------------
% Female susceptible
% ------------------------------------------------------------

SFout = uint16(65535 .* mat2gray(SfinalF));

imwrite( ...
    SFout, ...
    fullfile(outputFolder,'Susceptible_Females_Final.tif'));


% ------------------------------------------------------------
% Male susceptible
% ------------------------------------------------------------

SMout = uint16(65535 .* mat2gray(SfinalM));

imwrite( ...
    SMout, ...
    fullfile(outputFolder,'Susceptible_Males_Final.tif'));


% ------------------------------------------------------------
% Female infected
% ------------------------------------------------------------

IFout = uint16(65535 .* mat2gray(IfinalF));

imwrite( ...
    IFout, ...
    fullfile(outputFolder,'Infected_Females_Final.tif'));


% ------------------------------------------------------------
% Male infected
% ------------------------------------------------------------

IMout = uint16(65535 .* mat2gray(IfinalM));

imwrite( ...
    IMout, ...
    fullfile(outputFolder,'Infected_Males_Final.tif'));


%% ============================================================
% 23. FINAL SUMMARY
% ============================================================

fprintf('\n============================================\n');
fprintf(' OUTPUTS SAVED\n');
fprintf('============================================\n');

fprintf('\nOutput folder:\n');
fprintf('  %s\n',outputFolder);

fprintf('\nFinal populations:\n');

fprintf('  Susceptible females: %.4f\n', ...
    totalSF(end));

fprintf('  Susceptible males:   %.4f\n', ...
    totalSM(end));

fprintf('  Infected females:    %.4f\n', ...
    totalIF(end));

fprintf('  Infected males:      %.4f\n', ...
    totalIM(end));

fprintf('  Total susceptible:   %.4f\n', ...
    totalS(end));

fprintf('  Total infected:      %.4f\n', ...
    totalI(end));

fprintf('  Total deer:          %.4f\n', ...
    totalDeer(end));

fprintf('  Final prevalence:    %.4f%%\n', ...
    100 .* prevalence(end));

fprintf('\nFiles saved successfully.\n');


end
