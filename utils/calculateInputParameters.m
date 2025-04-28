function [projectileMassKg, caliberMm, baseAreaM2, initialPropellantMassKg, initialFreeVolumeM3] = ...
    calculateInputParameters(projectileWeightGr, projectileDiameterIn, propellantMassGr, propellantDensityKgm3, caseCapacityGrH2O, bulletLength_in, seatingDepth_in) % <- NUOVI INPUT

    % --- Conversion Constants ---
    KG_PER_GRAIN = 1 / 15432.35835;
    MM_PER_INCH = 25.4;
    METERS_PER_MM = 0.001;
    METERS_PER_INCH = 0.0254; % Nuovo
    WATER_DENSITY_KGM3 = 1000.0;
    M3_PER_GRAIN_H2O = KG_PER_GRAIN / WATER_DENSITY_KGM3;

    % --- Basic Input Validation ---
    % Aggiungere controlli per bulletLength_in, seatingDepth_in se necessario
    if projectileWeightGr <= 0 || projectileDiameterIn <= 0 || propellantMassGr < 0 || propellantDensityKgm3 <= 0 || caseCapacityGrH2O <= 0 || bulletLength_in <= 0 || seatingDepth_in <= 0
         error('calculateInputParameters: Input parameters must be positive (except potentially propellant mass).');
    end
    % ... (altri controlli)

    % --- Conversion and Derivation Calculations ---
    % 1. Projectile Mass
    projectileMassKg = projectileWeightGr * KG_PER_GRAIN;
    % 2. Projectile Caliber and Base Area
    caliberMm = projectileDiameterIn * MM_PER_INCH;
    projectileDiameterM = projectileDiameterIn * METERS_PER_INCH; % Nuovo: Diametro in metri
    projectileRadiusM = projectileDiameterM / 2.0;
    baseAreaM2 = pi * projectileRadiusM^2;
    % 3. Initial Propellant Mass
    initialPropellantMassKg = propellantMassGr * KG_PER_GRAIN;
    % 4. Initial Free Volume (V0) - MODIFICATO
    % 4a. Total internal case volume
    totalCaseVolumeM3 = caseCapacityGrH2O * M3_PER_GRAIN_H2O;
    % 4b. Volume occupied by solid propellant
    propellantSolidVolumeM3 = 0;
    if initialPropellantMassKg > 0
        propellantSolidVolumeM3 = initialPropellantMassKg / propellantDensityKgm3;
    end
    % 4c. Volume occupied by seated bullet (cylindrical approximation)
    seatingDepth_m = seatingDepth_in * METERS_PER_INCH;
    seatedBulletVolumeM3 = baseAreaM2 * seatingDepth_m;
    % 4d. Free volume V0
    initialFreeVolumeM3 = totalCaseVolumeM3 - propellantSolidVolumeM3 - seatedBulletVolumeM3;

    % --- Sanity Check for V0 ---
    if initialFreeVolumeM3 <= 1e-12 % Check if V0 is practically zero or negative
        warning('calculateInputParameters:VolumeNonPositive', ...
                ['Calculated initial volume V0 (%.3e m^3) is non-positive or very small.\n' ...
                 'Check case capacity (%.1f grH2O -> %.3e m^3), propellant charge (%.1f gr -> %.3e m^3), seating depth (%.3f in -> %.3e m^3 bullet vol), and density (%.0f kg/m^3).\n' ...
                 'Verify input data. Proceeding with V0 = %.3e m^3.'], ...
                initialFreeVolumeM3, caseCapacityGrH2O, totalCaseVolumeM3, propellantMassGr, propellantSolidVolumeM3, seatingDepth_in, seatedBulletVolumeM3, propellantDensityKgm3, initialFreeVolumeM3);
         % initialFreeVolumeM3 = max(initialFreeVolumeM3, 1e-9); % Opzionale: forza un minimo
    end

    % --- (Optional) Display Summary ---
    % ...
end