% Reset MATLAB
clear; close all;
FIGURE_COUNTER = 0;

% Define constants
planck_constant = 6.626 * 10 ^ -34; % units: J * s
electron_charge = 1.602 * 10 ^ -19; % units: C
electron_mass = 9.11 * 10 ^ -31; % units: kg
planck_eV = planck_constant / electron_charge; % units J * s * eV / J
speed_of_light = 3.0 * 10 ^ 8; % units: m / s
Boltzmann_constant = 1.381 * 10 ^ -23; % units: J / K
T = 300; % units: K

% Define main directory
Main_Directory = "C:\Users\fnfre\Desktop\Github\New Indoor PV/";

% Define material
Material_Name = "MoS_2";

% Define material's electrical parameters here
Material_BandGap = 1.27; % units: eV
Material_me = 0.71 * electron_mass;
Material_mh = 0.84 * electron_mass;

% Calculate effective conduction and valence band density of states and
% intrinsic carrier concentration
Material_NC = 10 ^ -6 * 2 * ((2 * pi * Material_me * Boltzmann_constant * T) / ...
    (planck_constant) ^ 2) ^ (3 / 2); % units: cm ^ -3
Material_NV = 10 ^ -6 * 2 * ((2 * pi * Material_mh * Boltzmann_constant * T) / ...
    (planck_constant) ^ 2) ^ (3 / 2); % units: cm ^ -3
Material_ni = sqrt(Material_NC * Material_NV) * exp(-(Material_BandGap * ...
    electron_charge) / (2 * Boltzmann_constant * T)); % units: cm ^ -3

% Define material's optical parameters
Material_Directory = strcat(Main_Directory, "Data/Material/", Material_Name, "/");
Material_k_File = strcat(Material_Directory, Material_Name, "-k.txt");
Material_k_Data = readmatrix(Material_k_File);
Material_n_File = strcat(Material_Directory, Material_Name, "-n.txt");
Material_n_Data = readmatrix(Material_n_File);

% Calculate alpha_2 and find n at the material's band gap energy
Material_Alpha2_Original = 4 * pi * Material_k_Data(:, 2:2) ./ Material_k_Data(:, 1:1); % units: nm ^ -1
Material_BandGap_Wavelength = (planck_eV * speed_of_light / ...
    Material_BandGap) * 10 ^ 9; % units: nm
Material_n = interp1(Material_n_Data(:, 1:1), Material_n_Data(:, 2:2), ...
    Material_BandGap_Wavelength); % no units
Material_4n2 = 4 * Material_n ^ 2; % 4 * n ^ 2, also no units
Material_n2 = Material_n ^ 2; % n ^ 2, also no units

% Define spectrum directory, light sources and intensities (in lux)
Spectrum_Directory = strcat(Main_Directory, "Data/Spectrum/");
Spectra_Names = ["LED"; "ESL"; "Halogen"; "Real AM15G"];
Intensities = [50; 150; 400; 500]; % units: lux

% Make sure that the material's alpha_2 data and spectrum data are as a function
% of the same wavelength domain
Minimum_Wavelength = 280.0;
Maximum_Wavelength = 4000.0;
New_Wavelength_Range = linspace(Minimum_Wavelength, Maximum_Wavelength, ...
    20 * (Maximum_Wavelength - Minimum_Wavelength) + 1).';
New_Energy_J_Range = 10 ^ 9 * planck_constant * speed_of_light./ ...
    New_Wavelength_Range; % units: (J * s * m / s) / (nm * 10 ^ -9 m / nm) = J / photon
New_Energy_eV_range = New_Energy_J_Range / electron_charge;
Material_Alpha2 = interp1(New_Wavelength_Range, Material_Alpha2_Original, ...
    New_Wavelength_Range); % no units

% Define recombination parameters
Material_Alpha1 = 0; % Free Carrier Absorption coefficient; units: nm ^ -1
Material_C = 10 ^ -29.66; % Auger coefficient; units: cm ^ 6 * s ^ -1
SRH_tau_range = [Inf 611e-9]; % units: s

% Supports up to eight tau_SRH values; if you would like to use more,
% add more colors to this matrix
% Colors = [[0.5 0 0]; [0.64 0.08 0.18]; [0.85 0.33 0.10];
%     [0.93 0.69 0.13]; [0 0.5 0]; [0 0.75 0.75];
%     [0.00 0.45 0.74]; [0.49 0.18 0.56]];
Colors_All = [[0.85 0.33 0.10]; [0.07 0.62 1.00];
    [0.25 0.50 0.07]; [0.49 0.18 0.56]];

% Define thickness range and specific thicknesses of interest to generate
% JV curves and spectral dependencies of the luminescent emission rate
L_Range = linspace(1, 1001, 1001).'; % units: nm
L_Interest = [100 1000]; % units: nm (put in ascending order)
for i = 1:length(L_Interest)
    assert(ismember(L_Interest(i), L_Range), ...
        'Your thickness of interest is not defined in your thickness range.');
end

% Time loggers to approximate time remaining
Time_Logger_Regular = 0;
Time_Logger_Interest = 0;
%% Operating Parameters

% Outer to Inner loop: Spectra --> Intensity --> Length --> tau_SRH
% Different data tables for each intensity. Each intensity folder has tables
% that are 1000 thicknesses * X tau_SRH values.

for spectra_index = 1:length(Spectra_Names)
    Spectra_Name = Spectra_Names(spectra_index);

    % Define figure output directories
    [~, ~, ~] = mkdir(strcat(Spectra_Name, "/Figures (.fig)"));
    Figures_fig_Directory = strcat(Spectra_Name, "/Figures (.fig)", "/");
    [~, ~, ~] = mkdir(strcat(Spectra_Name, "/Figures (.jpg)"));
    Figures_jpg_Directory = strcat(Spectra_Name, "/Figures (.jpg)", "/");
    [~, ~, ~] = mkdir(strcat(Spectra_Name, "/Figures by Intensity (.fig)"));
    Figures_Intensity_fig_Directory = strcat(Spectra_Name, "/Figures by Intensity (.fig)", "/");
    [~, ~, ~] = mkdir(strcat(Spectra_Name, "/Figures by Intensity (.jpg)"));
    Figures_Intensity_jpg_Directory = strcat(Spectra_Name, "/Figures by Intensity (.jpg)", "/");

    % Get the Intensity in W * m ^ -2
    Total_Intensities = readmatrix(strcat(Spectrum_Directory, Spectra_Name, ...
        "/", Spectra_Name, " Intensities and Luxes.txt"));

    for intensity_index = 1:length(Intensities)
        Intensity = Intensities(intensity_index); % units: Lux
        Total_Intensity = Total_Intensities(intensity_index, 1) * 10 ^ 3 * 10 ^ -4;
        % units: W * m ^ -2 = 10 ^ 3 * 10 ^ -4 * mW * cm ^ -2

        % Set up each spectra for each intensity, and make data output directories
        Spectrum_File = strcat(Spectrum_Directory, Spectra_Name, "/", Spectra_Name, " (", ...
            num2str(Intensity), " Lux).txt");
        Spectrum_Data = readmatrix(Spectrum_File);
        Spectrum = interp1(Spectrum_Data(:, 1:1), Spectrum_Data(:, 2:2), ...
            New_Wavelength_Range); % units: W * m ^ -2 * nm ^ -1
        clear Spectrum_Data;
        [~, ~, ~] = mkdir(strcat(Spectra_Name, "/", num2str(Intensity), ...
            " Lux Output Data"));
        Output_Data_Directory = strcat(Spectra_Name, "/", num2str(Intensity), ...
            " Lux Output Data/");

        % Section 1: Calculate efficiency limits with Shockley-Queisser (SQ) Model
        SQ_Voltage = linspace(0, 1.5, 1.5e3).'; % Voltage domain, units: V
        Black_Body_Photon_Flux = @(E_gap) (E_gap .^ 2) ./ (exp(E_gap ./ ((Boltzmann_constant ...
            / electron_charge) * T)) - 1); % units: E_gap in eV --> eV ^ 3
        
        % Radiative recombination
        RR0 = (10 ^ 3 * 10 ^ -4 * electron_charge * 2 * pi) / (speed_of_light ^ 2 ...
            * (planck_constant / electron_charge) ^ 3) * integral(Black_Body_Photon_Flux, ...
            Material_BandGap, Inf); % units: A / m ^ 2 = 10 ^ 3 * 10 ^ -4 * mA / cm ^ 2
        
        % Smallest wavelength (or largest energy) that will be absorbed
        Material_Lower_Bound_Wavelength = New_Wavelength_Range(1);
        % Largest wavelength (or smallest energy) that will be absorbed, up to the
        % band gap wavelength
        Material_Upper_Bound_Wavelength = Material_BandGap_Wavelength;
        % Wavelength range that will be absorbed
        Material_Wavelength_Domain = linspace(Material_Lower_Bound_Wavelength, ...
            Material_Upper_Bound_Wavelength, 20 * (Material_Upper_Bound_Wavelength ...
            - Material_Lower_Bound_Wavelength) + 1).';
        % Convert to wavelength to photon energy
        Material_Photon_Energy = 10 ^ 9 * planck_constant * speed_of_light ./ Material_Wavelength_Domain;
        % Select portion of spectrum that will be absorbed
        Material_Spectrum = interp1(New_Wavelength_Range, Spectrum, Material_Wavelength_Domain); % units: J / (nm * m ^ 2)
        
        % Calculate light-generated current by integrating the portion of the
        % spectrum that is absorbed
        SQ_L = 10 ^ 3 * 10 ^ -4 * electron_charge * trapz(Material_Wavelength_Domain, ...
            Material_Spectrum ./ Material_Photon_Energy); % (photon * C) / (nm * m ^ 2) = A / m ^ 2 =
            % 10 ^ 3 * 10 ^ -4 * mA / cm ^ 2
        
        % Find current density from radiative recombination and light-generated
        % current, using the ideal diode equation
        SQ_Current_Density = RR0 * (exp((electron_charge .* SQ_Voltage) ...
            ./ (Boltzmann_constant * T))) - SQ_L;
        
        % Calculate SQ parameters
        SQ_Jsc = abs(interp1(SQ_Voltage, SQ_Current_Density, 0));
        SQ_Current_Density_new = SQ_Current_Density;
        [SQ_Current_Density_new, index] = unique(SQ_Current_Density_new);
        SQ_Voc = interp1(SQ_Current_Density_new, SQ_Voltage(index), 0); % V
        SQ_FF = abs(min(SQ_Current_Density .* SQ_Voltage)) / abs(SQ_Voc * SQ_Jsc);
        SQ_Pmax = abs(SQ_Jsc * SQ_Voc * SQ_FF); % mW / cm ^ 2
        SQ_eff = SQ_Pmax / Total_Intensity * 100;
        
        % Export the SQ data
        writematrix([SQ_Voltage SQ_Current_Density], strcat(Output_Data_Directory, ...
            "Shockley-Queisser Model JV Curve.txt"), 'Delimiter', 'tab');
        writematrix([SQ_Voc; SQ_Jsc; SQ_FF; SQ_Pmax; SQ_eff], strcat(Output_Data_Directory, ...
            "Shockley-Queisser Parameters.txt"), 'Delimiter', 'tab');

        % Section 2: Calculate efficiency limits with extended Tiedje-Yablonovitch
        % (TY) Model (includes Auger recombination) and defect-assisted
        % Shockley-Reed-Hall (SRH) recombination
        
        % Data output tables
        % Parameters
        Voc = zeros(length(L_Range), length(SRH_tau_range));
        Jsc = zeros(length(L_Range), length(SRH_tau_range));
        FF = zeros(length(L_Range), length(SRH_tau_range));
        Vmax = zeros(length(L_Range), length(SRH_tau_range));
        Imax = zeros(length(L_Range), length(SRH_tau_range));
        Pout = zeros(length(L_Range), length(SRH_tau_range));
        eff = zeros(length(L_Range), length(SRH_tau_range));
        
        % Recombination components
        AugerComponent = zeros(length(L_Range), length(SRH_tau_range));
        FreeCarrierComponent = zeros(length(L_Range), length(SRH_tau_range));
        InternalLuminescenceComponent = zeros(length(L_Range), length(SRH_tau_range));
        ReabsorptionComponent = zeros(length(L_Range), length(SRH_tau_range));
        ExternalEmissionComponent = zeros(length(L_Range), length(SRH_tau_range));
        SRHComponent = zeros(length(L_Range), length(SRH_tau_range));
        
        % Lifetimes
        Auger_lifetime = zeros(length(L_Range), length(SRH_tau_range));
        Radiative_lifetime = zeros(length(L_Range), length(SRH_tau_range));
        SRH_lifetime = zeros(length(L_Range), length(SRH_tau_range));

        % Supplementary Note 1, Equation 3
        Black_Body_Photon_Flux = (2 * Material_n2) / (speed_of_light ^ 2 * ...
            planck_eV ^ 3) * (New_Energy_eV_range .^ 2) ./ (exp(New_Energy_eV_range ...
            ./ ((Boltzmann_constant / electron_charge) * T)) - 1);
        
        thickness_interest_index = 0;
        for thickness_index = 1:length(L_Range)
        
            % Print current spectrum, intensity, thickness, and time
            disp(strcat("Current Spectrum: ", Spectra_Name, ", ", num2str(Intensity), " Lux"));
            L = L_Range(thickness_index);
            fprintf(['Current Thickness: ' num2str(L), ' nm (Started), ']);
            Time1 = datetime('now');
            fprintf(['Current Time Now: ' datestr(Time1) '.\n']);
        
            if ismember(L, L_Interest)
                [~, ~, ~] = mkdir(strcat(Output_Data_Directory, "L = ", num2str(L), " nm Data"));
                L_Interest_Directory = strcat(Output_Data_Directory, "L = ", num2str(L), " nm Data/");
            end
        
            % Supplementary Note 1, Equation 2
            Absorbance = Material_Alpha2 ./ (Material_Alpha2 + Material_Alpha1 + (1 ./ (Material_4n2 .* L)));
            
            f_range = [0; linspace(0.8, 1, 100).'];
            if ismember(L, L_Interest) % use whole range for thicknesses of interest
                f_range = linspace(0, 1, 1000).';
            end
            
            Current_Density = zeros(length(f_range), length(SRH_tau_range));
            Voltage = zeros(length(f_range), length(SRH_tau_range));
            
            for tau_index = 1:length(SRH_tau_range)
                SRH_tau = SRH_tau_range(tau_index);
        
                % Supplementary Note 1, Equation 4
                Jsc(thickness_index, tau_index) = 10 ^ 3 * 10 ^ -4 * electron_charge *...
                    trapz(New_Wavelength_Range, Spectrum .* Absorbance ./ New_Energy_J_Range);
                % integrand units: (J / s) / (nm * m ^ 2) / J = 1 / (nm * m ^ 2 * s)
                % units: (photon * C) / (nm * m ^ 2 * s) = A / m ^ 2 = 10 ^ 3 * 10 ^ -4 * mA / cm ^ 2
        
                for fraction_index = 1:length(f_range)
                    f = f_range(fraction_index);
        
                    syms x;
                    % Supplementary Note 1, Equation 11
                    % units (first term): (nm ^ -1 + nm ^ -1) * 10 ^ 7 nm / cm = 1 / cm
                    % units (second term): (10 ^ -4 m ^ 2 / cm ^ 2) * (1 / (s * m ^ 2 * eV)) * eV = cm ^ -2 / s
                    % units (third term): cm ^ 6 * s ^ -1 * cm ^ -9 = cm ^ -3 / s
                    % units (right hand side): ((C / s) / cm ^ 2) / (C * nm * 10 ^ -7 cm / nm) = cm ^ -3 / s
                    if SRH_tau == Inf % No SRH
                        eqn_RadiativeAugerFreeSRH = 4 * pi * (10 ^ 7 * (Material_Alpha1 * ...
                            exp(x / 2) + 1 / (Material_4n2 .* L))) * exp(x) * 10 ^ -4 * ...
                            trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* ...
                            Absorbance)) + Material_C * Material_ni ^ 3 * exp((3 / 2) * ...
                            x) == (Jsc(thickness_index, tau_index) * 10 ^ -3) / ...
                            (electron_charge * L * 10 ^ -7) * (1 - f);
                    else
                        eqn_RadiativeAugerFreeSRH = 4 * pi * (10 ^ 7 * (Material_Alpha1 * ...
                            exp(x / 2) + 1 / (Material_4n2 .* L))) * exp(x) * 10 ^ -4 * ...
                            trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* ...
                            Absorbance)) + Material_C * Material_ni ^ 3 * exp((3 / 2) * ...
                            x) + Material_ni * exp (x / 2) / SRH_tau == ...
                            (Jsc(thickness_index, tau_index) * 10 ^ -3) / (electron_charge * ...
                            L * 10 ^ -7) * (1 - f);
                    end
        
                    X_RadiativeAugerFreeSRH = vpasolve(eqn_RadiativeAugerFreeSRH, x);
        
                    Current_Density(fraction_index, tau_index) = f * Jsc(thickness_index, tau_index);
                    if f == 1
                        if any(Voltage(:, tau_index) < 0)
                            Voltage(fraction_index, tau_index) = Voltage(fraction_index - 1, tau_index);
                            Jsc(thickness_index, tau_index) = interp1(Voltage(1:end - 1, tau_index), ...
                                Current_Density(1:end - 1, tau_index), 0);
                        else
                            Voltage(fraction_index, tau_index) = 0;
                        end
                    else
                        Voltage(fraction_index, tau_index) = (Boltzmann_constant * ...
                            T * X_RadiativeAugerFreeSRH) / (electron_charge);
                    end
        
                    if f == 0
                        Voc(thickness_index, tau_index) = (Boltzmann_constant * T * ...
                            X_RadiativeAugerFreeSRH) / (electron_charge);
                    end
                end

                if ismember(L, L_Interest)
                    writematrix([Voltage Current_Density], strcat(L_Interest_Directory, ...
                        "JV Data at ", num2str(L), " nm.txt"), 'Delimiter', 'tab');
                end
        
                % Update the table with parameters
                FF(thickness_index, tau_index) = max(Current_Density(:, tau_index) .* Voltage(:, tau_index)) /...
                    (Voc(thickness_index, tau_index) * Jsc(thickness_index, tau_index));
                for max_index = 1:length(Voltage(:, tau_index))
                   if Voltage(max_index, tau_index) * Current_Density(max_index, tau_index) == max(Voltage(:, tau_index) .* Current_Density(:, tau_index))
                       Vmax(thickness_index, tau_index) = Voltage(max_index, tau_index);
                       Imax(thickness_index, tau_index) = Current_Density(max_index, tau_index);
                       Pout(thickness_index, tau_index) = Voltage(max_index, tau_index) * Current_Density(max_index, tau_index);
                   end
                end
                eff(thickness_index, tau_index) = max(Current_Density(:, tau_index) .* Voltage(:, tau_index)) / Total_Intensity * 100;
        
                % units: mA / cm ^ 2
                AugerComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * Material_C * Material_ni ^ 3 * ...
                    exp((3 / 2) * ((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)));
                FreeCarrierComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * 4 * pi * 10 ^ -4 * ...
                    exp((3 / 2) * ((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T))) * ...
                    (trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 * Material_Alpha1)));
                InternalLuminescenceComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * 4 * pi * 10 ^ -4 * ...
                    exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)) * ...
                    (Material_Alpha1 * exp((1 / 2) * (electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)) * ...
                    trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7)) + ...
                    trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 .* Material_Alpha2)) + ...
                    (1 / (Material_4n2 .* L * 10 ^ -7)) * trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance)));
                ReabsorptionComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * 4 * pi * 10 ^ -4 * ...
                    exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)) * ...
                    (trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 .* Material_Alpha2)));
                ExternalEmissionComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * 4 * pi * 10 ^ -4 * ...
                    exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)) * ...
                    (1 / (Material_4n2 .* L * 10 ^ -7)) * trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance));
                SRHComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * Material_ni * ...
                    exp ((1 / 2) * ((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T))) / SRH_tau;
               
                % units: s
                Auger_lifetime(thickness_index, tau_index) = 1 / (Material_C * Material_ni ^ 2 * exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)));
                Radiative_lifetime(thickness_index, tau_index) = 1 ./ ((Material_ni * sqrt(exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)))) .* ...
                    (((4 * pi * 10 ^ -4 * (trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 * Material_Alpha1)) + ...
                    trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 .* Material_Alpha2)) + ...
                    (1 / (Material_4n2 .* L * 10 ^ -7)) * trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance))))) ./ ...
                    (Material_ni ^ 2)));
                SRH_lifetime(thickness_index, tau_index) = SRH_tau;
            end
        
            % Print current thickness and time
            fprintf(['Current Thickness: ' num2str(L), ' nm (Finished), ']);
            Time2 = datetime('now');
            fprintf(['Current Time Now: ' datestr(Time2) '.\n']);
        
            if ismember(L, L_Interest)
                Time_Logger_Interest = Time2 - Time1;
        
                % Move to the next thickness of interest
                thickness_interest_index = thickness_interest_index + 1;
            else
                Time_Logger_Regular = Time2 - Time1;
            end
        
            % Print number of thicknesses remaining
            remaining_lengths = length(L_Range) - thickness_index;
            remaining_lengths_interest = length(L_Interest) - thickness_interest_index;
            fprintf(['Thicknesses Remaining: ' num2str(remaining_lengths) ', ']);

            % Print number of thicknesses remaining
            remaining_spectra = length(Spectra_Names) - spectra_index;
            fprintf(['Spectra Remaining: ' num2str(remaining_spectra) ', ']);
            
            % Print number of thicknesses remaining
            remaining_intensities = length(Intensities) - intensity_index;
            fprintf(['Intensities Remaining: ' num2str(remaining_intensities) '.\n']);
            
            % Print approximate time remaining
            if thickness_interest_index > 0
                % More accurate prediction based on time it takes to calculate
                % thicknesses of interests
                remaining_time = (remaining_spectra * length(Intensities) * (length(L_Range) - length(L_Interest)) + ...
                    remaining_intensities * (length(L_Range) - length(L_Interest)) + ...
                    remaining_lengths - remaining_lengths_interest) * Time_Logger_Regular + ...
                    (remaining_spectra * length(Intensities) * (length(L_Interest)) + ...
                    remaining_intensities * (length(L_Interest)) + ...
                    remaining_lengths_interest) * Time_Logger_Interest;

            else
                remaining_time = (remaining_spectra * length(Intensities) * length(L_Range) + ...
                    remaining_intensities * length(L_Range) + ...
                    remaining_lengths) * Time_Logger_Regular;
            end
            [hours, minutes, seconds] = hms(remaining_time);
            fprintf(['Approximate Time Remaining: ' num2str(hours) ' Hours, ' ...
                num2str(minutes) ' Minutes, ' num2str(seconds) ' Seconds\n']);
        end

        % Export the data
        writematrix([L_Range Voc], strcat(Output_Data_Directory, "Voc.txt"), 'Delimiter', 'tab');
        writematrix([L_Range Jsc], strcat(Output_Data_Directory, "Jsc.txt"), 'Delimiter', 'tab');
        writematrix([L_Range FF], strcat(Output_Data_Directory, "FF.txt"), 'Delimiter', 'tab');
        writematrix([L_Range Vmax], strcat(Output_Data_Directory, "Vmax.txt"), 'Delimiter', 'tab');
        writematrix([L_Range Imax], strcat(Output_Data_Directory, "Imax.txt"), 'Delimiter', 'tab');
        writematrix([L_Range Pout], strcat(Output_Data_Directory, "Pout.txt"), 'Delimiter', 'tab');
        writematrix([L_Range eff], strcat(Output_Data_Directory, "eff.txt"), 'Delimiter', 'tab');
        
        writematrix([L_Range AugerComponent], strcat(Output_Data_Directory, "Component_Auger.txt"), 'Delimiter', 'tab');
        writematrix([L_Range FreeCarrierComponent], strcat(Output_Data_Directory, "Component_FreeCarrier.txt"), 'Delimiter', 'tab');
        writematrix([L_Range InternalLuminescenceComponent], strcat(Output_Data_Directory, "Component_InternalLuminescence.txt"), 'Delimiter', 'tab');
        writematrix([L_Range ReabsorptionComponent], strcat(Output_Data_Directory, "Component_Reabsorption.txt"), 'Delimiter', 'tab');
        writematrix([L_Range ExternalEmissionComponent], strcat(Output_Data_Directory, "Component_ExternalEmission.txt"), 'Delimiter', 'tab');
        writematrix([L_Range SRHComponent], strcat(Output_Data_Directory, "Component_SRH.txt"), 'Delimiter', 'tab');
        
        writematrix([L_Range Auger_lifetime], strcat(Output_Data_Directory, "Lifetime_Auger.txt"), 'Delimiter', 'tab');
        writematrix([L_Range Radiative_lifetime], strcat(Output_Data_Directory, "Lifetime_Radiative.txt"), 'Delimiter', 'tab');
        writematrix([L_Range SRH_lifetime], strcat(Output_Data_Directory, "Lifetime_SRH.txt"), 'Delimiter', 'tab');
    
        % Section 3: Plotting the data (Open Circuit Voltage, Short Circuit
        % Current Density, Fill Factor, Output Power, and Power Conversion
        % Efficiency) for each intensity individually
        
        SQ_Parameters = readmatrix(strcat(Output_Data_Directory, ...
            "Shockley-Queisser Parameters.txt"));
        
        Text_Filename = ["Voc"; "Jsc"; "FF"; "Pout"; "eff"];
        
        Y_Labels = ["V_{oc} (V)"; "J_{sc} (mA cm^{-2})"; ...
            "FF (%)"; "P_{out} (mW cm^{-2})"; "PCE (%)"];
        
        Titles = ["Open Circuit Voltage"; "Short Circuit Current Density"; ...
            "Fill Factor"; "Output Power"; "Power Conversion Efficiency"];
        
        for i = 1:length(Text_Filename)
            Parameters = readmatrix(strcat(Output_Data_Directory, Text_Filename(i), ".txt"));
            FIGURE_COUNTER = FIGURE_COUNTER + 1;
            figure(FIGURE_COUNTER);
        
            hold on;
            for j = 1:length(SRH_tau_range)
                semilogx(Parameters(:, 1), Parameters(:, j + 1), '-', 'color', Colors_All(j, :), 'linewidth', 3);
            end
            semilogx([Parameters(1, 1) Parameters(end, 1)], [SQ_Parameters(i) SQ_Parameters(i)], '-', 'color', '[0 0 0]', 'linewidth', 3);
        
            xlabel('Thickness (nm)');
            ylabel(Y_Labels(i));
            hTitle = title(strcat(Material_Name, " ", Titles(i), " vs. Thickness"));
            set(gca, 'FontSize', 30);
            set(hTitle, 'FontSize', 30);
            set(gca, 'XScale', 'log');
            set(gca, 'linewidth', 2);
            xlim([5 1000]);
            
            set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
            pbaspect([1.2 1 1]);
            set(gcf, 'PaperOrientation', 'landscape');
            set(gcf, 'PaperUnits', 'normalized');
            set(gcf, 'PaperPosition', [0 0 1 1]);
            box on;
        
            savefig(strcat(Figures_Intensity_fig_Directory, num2str(Intensity), " Lux ", Titles(i), " vs. Thickness.fig"));
            saveas(gcf, strcat(Figures_Intensity_jpg_Directory, num2str(Intensity), " Lux ", Titles(i), " vs. Thickness.jpg"));
            close;
        end
    end

    % Section 4: Plotting the data (Open Circuit Voltage, Short Circuit
    % Current Density, Fill Factor, Output Power, and Power Conversion
    % Efficiency) and J-V curves for all intensities

    for a = 1:length(Text_Filename)
        FIGURE_COUNTER = FIGURE_COUNTER + 1;
        figure(FIGURE_COUNTER);
        hold on;

        for b = 1:length(Intensities)
            Intensity = Intensities(b);

            Data_Directory = strcat(Spectra_Name, "/", num2str(Intensity), ...
                " Lux Output Data/");

            SQ_Parameters = readmatrix(strcat(Data_Directory, ...
                "Shockley-Queisser Parameters.txt"));

            Parameters = readmatrix(strcat(Data_Directory, Text_Filename(a), ".txt"));

            for c = 1:length(SRH_tau_range)
                if mod(c, 2) == 1
                    LineType = ':';
                else
                    LineType = '-.';
                end

                semilogx(Parameters(:, 1), Parameters(:, c + 1), LineType, 'color', Colors_All(b, :), 'linewidth', 3);
                hold on;
            end
            semilogx([Parameters(1, 1) Parameters(end, 1)], [SQ_Parameters(a) SQ_Parameters(a)], '-', 'color', Colors_All(b, :), 'linewidth', 3);
        end
        xlabel('Thickness (nm)');
        ylabel(Y_Labels(a));
        hTitle = title(strcat(Material_Name, " ", Titles(a), " vs. Thickness"));
        set(gca, 'FontSize', 30);
        set(hTitle, 'FontSize', 30);
        set(gca, 'linewidth', 2);
        set(gca, 'XScale', 'log');
        xlim([5 1000]);
    
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
        pbaspect([1.2 1 1]);
        set(gcf, 'PaperOrientation', 'landscape');
        set(gcf, 'PaperUnits', 'normalized');
        set(gcf, 'PaperPosition', [0 0 1 1]);
        box on;


        savefig(strcat(Figures_fig_Directory, Titles(a), " vs. Thickness.fig"));
        saveas(gcf, strcat(Figures_jpg_Directory, Titles(a), " vs. Thickness.jpg"));
        close;
    end

    % Section 5: Plotting the JV curve for all intensities
    for d = 1:length(L_Interest)
        Length = L_Interest(d);

        FIGURE_COUNTER = FIGURE_COUNTER + 1;
        figure(FIGURE_COUNTER);

        for e = 1:length(Intensities)
            Intensity = Intensities(e);

            Data_Directory = strcat(Spectra_Name, "/", num2str(Intensity), ...
                " Lux Output Data/");

            SQ_JV = readmatrix(strcat(Data_Directory, ...
                    "Shockley-Queisser Model JV Curve.txt"));
            SQ_Parameters = readmatrix(strcat(Data_Directory, ...
                "Shockley-Queisser Parameters.txt"));

            JV_Curve = readmatrix(strcat(Data_Directory, "L = ", num2str(Length), ...
                " nm Data/JV Data at ", num2str(Length), " nm.txt"));
    
            for g = 1:length(SRH_tau_range)
                if mod(g, 2) == 1
                    LineType = ':';
                else
                    LineType = '-.';
                end
                plot(JV_Curve(:, g), JV_Curve(:, length(SRH_tau_range) + g), LineType, 'color', Colors_All(e, :), 'linewidth', 3);
                hold on;
            end
            plot(SQ_JV(:, 1), -SQ_JV(:, 2), '-', 'color', Colors_All(e, :), 'linewidth', 3);
        end
        xlabel('Voltage (V)');
        ylabel('Current density (mA cm^{-2})');
        hTitle = title(strcat("JV Curve, Thickness = ", num2str(Length), " nm"));
        set(gca, 'FontSize', 30);
        set(hTitle, 'FontSize', 30);
        set(gca, 'linewidth', 2);
        xlim([-0.1 * SQ_Parameters(1) 1.25 * SQ_Parameters(1)]);
        ylim([-0.1 * SQ_Parameters(2) 1.25 * SQ_Parameters(2)]);
        
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
        pbaspect([1.2 1 1]);
        set(gcf, 'PaperOrientation', 'landscape');
        set(gcf, 'PaperUnits', 'normalized');
        set(gcf, 'PaperPosition', [0 0 1 1]);
        box on;

        % Export the plot and the data
        savefig(strcat(Figures_fig_Directory, "JV Curve, ", num2str(Length), " nm.fig"));
        saveas(gcf, strcat(Figures_jpg_Directory, "JV Curve, ", num2str(Length), " nm.jpg"));
        close;
    end
end
