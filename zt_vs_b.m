function zt_vs_b()
%zt_vs_b Plots ZT vs. b and locate measurement data on the plot.
%
%   This function accepts measured thermoelectric data and places it on the
%   maximum ZT vs. b map. The load line for the measured data against 
%   maximum ZT vs. b is also shown. 
%   -----------------------------------------------------------------------
%   Inputs
%
%   The program takes measured thermoelectric quantities being either
%   1. Seebeck, Resistivity, and ZT
%   2. Seebeck, Resistivity, and Total Thermal Conductivity
%   3. Seebeck, Resistivity, and Lattice Thermal Conductivity
%   
%   Notice that ZT and Thermal conductivity cannot be given to the program at the same time. 
%   The program will use one and calculate the other.
%
%   For details, see "Input parameters from experimental measurements" section below
%   -----------------------------------------------------------------------
%   Method synopsis
%
%   1. Using the given Seebeck coefficient, the Fermi level is determined from a single parabolic
%      band model. 
%   2. Using the given resistivity and the determined Fermi level, the resistivity and electronic 
%      thermal conductivity at different Fermi energies are calculated.
%   3. Using the given ZT or thermal conductivity, all the rest thermoelectric quantities of interest
%      can be calculated.
%   4. Vary the lattice thermal conductivity and locate the optimal Fermi energy that maximize the ZT
%      for each lattice thermal conductivity value. Covert this into a ZT|max vs. b curve.
%   -----------------------------------------------------------------------
%   Outputs
%
%   Figure #1: ZT|_{max} vs. b_l
%   Figure #2: ZT|_{max} vs. b_total
%   Figure #3: Seebeck coefficient vs. b
%   Figure #4: Electrical resistivity vs. b
%   Figure #5: Thermal conductivity vs. b
%   Figure #6: Lorenz number vs. b
%   -----------------------------------------------------------------------
%   This script follows the general procedure published in
%   "Transport property analysis method for thermoelectric materials: material quality factor and the effective mass model"
%   SD Kang, GJ Snyder - arXiv preprint arXiv:1710.06896, 2017
%
%   This function is written by
%   Authors: Xufeng Wang, Evan Witkosk, and Mark Lundstrom
%   Contact: lundstro@ecn.purdue.edu
%   Copyright 2017 Purdue University

close all; clear; clc

%% Input parameters from experimental measurements
% Required inputs
T             = 300; % (K) Temperature

seebeck_coeff = 100e-6; % (V/K) Seebeck coefficient
resistivity   = 1e-3; % (Ohm-cm) Electrical resistivity

% Either specify ZT or thermal conductivity, NEVER both
ZT            = []; % Thermoelectric Figure of Merit

kappa_option  = 'lattice'; % ('lattice'/'total') use either lattice thermal conductivity or total thermal conductivity
kappa         = 0.5; % (W/K-m) Thermal conductivity (total or lattice)

%% Simulation parameters
% Numerical parameters
b_range = [1e-3 0.4]; % desired range for b_L
b_desired_number = 100; % desired number of data points for b within range

tol = 1e-7; % integration boundary for the Fermi window. Integration terminates when integrant value reach peak*tol.

energy_res = 1e-4; % (eV) energy grid resolution for integration
fermi_list = -1:energy_res:1; % position of the fermi level with respect to band edge. E_C = 0 eV

r_coeff = 0; % () energy dependence of the Mean Free Path(mfp). mfp(E) = mfp_o*(E-E_C)^r_coeff

% Physical constants
q  = 1.6e-19      ; % (C) elementary charge
h  = 4.135667e-15 ; % (eV-s) Planck constant
kb = 8.617e-5     ; % (eV/K) Boltzmann constant

%% Preparations and pre-check
if isempty(ZT) && isempty(kappa) || ~isempty(ZT) && ~isempty(kappa)
    error('Must specifiy either the ZT or the thermal conductivity, but never both.');
end

seebeck_coeff = abs(seebeck_coeff); % Take the absolute value of Seebeck coefficient.

if resistivity <= 0; error('Resistivity must be more than zero'); end
if T <= 0; error('T must be more than zero'); end
if ~isempty(kappa); if kappa <= 0; error('kappa must be more than zero'); end; end
if ~isempty(ZT); if ZT <= 0; error('ZT must be more than zero'); end; end
 
%% Comput thermoelectric properties
fprintf('Computing thermoelectric properties ...\n');

% Preallocation
lorenz_list  = zeros(1, length(fermi_list));
seebeck_list = zeros(1, length(fermi_list));
e_cond_list  = zeros(1, length(fermi_list));

% find Fermi window integration range
fun = @(x) exp(x./(kb.*T))./(1+exp(x./(kb.*T))).^2 - tol;
f_window_range = abs(fzero(fun, 0)); % This f_window_range gives a Fermi tail less than tol

% normalized TDF
tdf_function(:,1) = (fermi_list(1)-f_window_range):energy_res:(fermi_list(end)+f_window_range);
tdf_function(tdf_function(:,1)>=0,2) = tdf_function(tdf_function(:,1)>=0,1).^(r_coeff+1); % main
tdf_function(tdf_function(:,2) == Inf,2) = 0; % remove infinity from r_coeff < 0 cases near band edge

for index = 1:length(fermi_list)
    % progress bar
    if rem(index,round(length(fermi_list)/10)) == 0
        fprintf('..%d%%',10*floor(index/round(length(fermi_list)/10)));
    end
    
    % find energy range for integration
    energy_lower_bnd = max([min(tdf_function(:,1)) ...
        max([0 fermi_list(index)-f_window_range])]); % Fermi energy grid upper bound initial guess
    energy_upper_bnd = min([max(tdf_function(:,1)) ...
        max([f_window_range fermi_list(index)+f_window_range])]); % Fermi energy grid upper bound initial guess
    
    % ensure energy grid is fine enough
    if energy_res < (energy_upper_bnd-energy_lower_bnd)/1000
        % if specified energy_res gives MORE than 1000 grid point within
        % Fermi window, use energy_res
        energy_list = energy_lower_bnd:energy_res:energy_upper_bnd;
    else
        % if specified energy_res gives LESS than 1000 grid points within
        % Fermi window, generate 1000 grid points instead
        energy_list = linspace(energy_lower_bnd, energy_upper_bnd, 1000);
    end
    
    d_fermi = (1/(kb*T)).*exp((energy_list-fermi_list(index))./(kb.*T))./(1+exp((energy_list-fermi_list(index))./(kb.*T))).^2; % Fermi window
    
    this_tdf_function = interp1(tdf_function(:,1), tdf_function(:,2), energy_list); % normalized transport distribution function
    
    e_cond_resolved = (2*q/h).*d_fermi.*this_tdf_function; % differential conductivity
    
    e_cond = trapz(energy_list, e_cond_resolved, 2);
    e_cond_list(index) = e_cond; % (arb.) conductivity
    
    e_cond_e = trapz(energy_list, e_cond_resolved.*-1.*(energy_list-fermi_list(index)), 2);
    e_cond_e2 = trapz(energy_list, e_cond_resolved.*(energy_list-fermi_list(index)).^2, 2);
    
    lorenz_list(index) = e_cond_e2./e_cond - (e_cond_e./e_cond).^2;
    lorenz_list(index) = lorenz_list(index)/(kb*T)^2; % Normalized Lorenz number
    
    seebeck_list(index) = abs(e_cond_e./e_cond./T); % (V/K) Seebeck coefficient
end

fprintf('..Done!\n');

% remove nan and Inf values
% this occurs most likely because Fermi level being too far from band edge
removed_index = ~isfinite(seebeck_list);
fermi_list(removed_index) = [];
lorenz_list(removed_index) = [];
seebeck_list(removed_index) = [];
e_cond_list(removed_index) = [];

% check if measured Seebeck coefficient is within the calculated range
if abs(seebeck_coeff) < min(seebeck_list) || abs(seebeck_coeff) > max(seebeck_list)
    error('The input measured Seebeck coefficient is outside of calculated by given Fermi energy range. Modify the Fermi energy range to correct this error.');
end

%% Locate quantities based on measurement values
% Fermi level at measured Seebeck coefficient
fermi_at_seebeck = interp1(seebeck_list, fermi_list, abs(seebeck_coeff));

% scale conductivity vs. Fermi level, so it gives same conductivity as
% measured at the same Fermi level
e_cond_at_seebeck = interp1(fermi_list, e_cond_list, fermi_at_seebeck);
e_cond_list = 1/resistivity/e_cond_at_seebeck.*e_cond_list.*1e2; % (arb. -> S/m)

% calculate electronic thermal conductivity
e_kappa_list = lorenz_list.*(kb)^2.*e_cond_list.*T; % (W/K-m) Electronic thermal conducitvity
e_kappa_at_seebeck = interp1(fermi_list, e_kappa_list, fermi_at_seebeck); % (W/K-m) predicted measured electrical thermal conductivity

% calculate the measured ZT
if isempty(ZT)
    % kappa is given as input
    
    switch kappa_option
        case 'total'
            % input kappa is total thermal conductivity
            target_kappa_l = kappa - e_kappa_at_seebeck;
            
        case 'lattice'
            % input kappa is lattice thermal conductivity
            target_kappa_l = kappa;
            kappa = target_kappa_l + e_kappa_at_seebeck; % rename "kappa" as the total thermal conductivity
            
        otherwise
            error('Invalid thermal conductivity option. Can only be "total" or "lattice"')
    end
    
    ZT = 1/resistivity*1e2*seebeck_coeff.^2*T./kappa;
    
else
    % ZT is given as input
    kappa = 1/resistivity*1e2*seebeck_coeff.^2*T./ZT;
    target_kappa_l = kappa - e_kappa_at_seebeck;
end

if kappa <= e_kappa_at_seebeck
    % input kappa is total thermal conductivity and it is less than kappa_e
    error('Calculated electrical thermal conductivity %1.2f W/K-m is higher than the input total %1.2f W/K-m.', e_kappa_at_seebeck, kappa);
end

fprintf('--------------------------------------------------------\n');
fprintf('Input material parameters summary\n\n')
fprintf('ZT = %3f at T = %3f K\n\n', ZT, T);
fprintf('Seebeck Coefficient             = %1.2e uV/K;  Resistivity                  = %1.2e Ohm-cm\n\n', seebeck_coeff*1e6, resistivity);
fprintf('Total thermal conductivity      = %1.2e W/K-m;\n', target_kappa_l+e_kappa_at_seebeck);
fprintf('Electronic thermal conductivity = %1.2e W/K-m; Lattice thermal conductivity = %1.2e W/K-m\n', e_kappa_at_seebeck, target_kappa_l);
fprintf('Lorenz number                   = %1.2f\n', interp1(fermi_list, lorenz_list, fermi_at_seebeck));
fprintf('--------------------------------------------------------\n');

% the range for b should be no less than specified b_min and (experimentally measure b * 0.5) and no more than b_max and (experimentally measured b * 1.5)
target_b = T*kb^2/resistivity/target_kappa_l*1e2; % input b_l
target_b_tot = ZT/(seebeck_coeff/kb)^2; % input b_total

if target_b*1.5 > max(b_range)
    b_range = [min(b_range) target_b*1.5];
    fprintf('Notice: b range reset to %e to %e to cover input values.\n', min(b_range), max(b_range))
end
if target_b*0.5 < min(b_range)
    b_range = [target_b*0.5 max(target_b)];
    fprintf('Notice: b range reset to %e to %e to cover input values.\n', min(b_range), max(b_range))
end
if min(b_range) < 1e-8
    % avoid division by zero problem or b erroursly set to below zero
    b_range = 1e-8;
    fprintf('Notice: b range reset to %e to %e to avoid zero/negative values.\n', min(b_range), max(b_range))
end
if min(b_range) == max(b_range)
    error('b range minimum and maximum values cannot be the same.')
end

fprintf('Computing ZT vs. b ...\n');

%% Compute ZT vs. b by varying lattice conductivity

b_within_range = 0; % number of b data points within desired b range

kappa_l_lower = 1./max(b_range).*T.*kb^2.*min(e_cond_list); % initial guess for a kappa_l that meets b lower limit
kappa_l_upper = 1./min(b_range).*T.*kb^2.*max(e_cond_list); % initial guess for a kappa_l that meets b upper limit

while b_within_range < b_desired_number
    % loop until desired number of b points are found
    kappa_l_list = logspace(log10(kappa_l_lower),log10(kappa_l_upper), b_desired_number*2); % kappa_l
    
    % Preallocation
    b_resolved       = zeros(length(kappa_l_list), length(fermi_list)); % Fermi energy resolved b_L at each kappa_l
    b_tot_resolved   = zeros(length(kappa_l_list), length(fermi_list)); % Fermi energy resolved b_total at each kappa_l
    zt_resolved      = zeros(length(kappa_l_list), length(fermi_list)); % Fermi energy resolved ZT at each kappa_l
    
    zt_at_max     = zeros(1,length(kappa_l_list)); % max ZT at each kappa_l
    b_at_max      = zeros(1,length(kappa_l_list)); % b_L at max ZT at each kappa_l
    b_tot_at_max  = zeros(1,length(kappa_l_list)); % b_total at max ZT at each kappa_l 
    s_at_max      = zeros(1,length(kappa_l_list)); % Seebeck coefficient at max ZT at each kappa_l
    e_cond_at_max = zeros(1,length(kappa_l_list)); % conductivity at max ZT at each kappa_l
    l_at_max      = zeros(1,length(kappa_l_list)); % Lorenz number at max ZT at each kappa_l
    
    for kappa_l_index = 1:length(kappa_l_list)
        % for each kappa_l
        this_kappa_l = kappa_l_list(kappa_l_index);
        this_b = T*kb^2.*e_cond_list./this_kappa_l;
        this_zt = (seebeck_list./kb).^2./(lorenz_list+1./this_b);
        this_b_tot = this_zt./(seebeck_list./kb).^2;
        
        [~, index_at_max_zt] = max(this_zt); % locate a rough maximum ZT
        
        if index_at_max_zt == 1 || index_at_max_zt == length(this_zt)
            % data is no longer reliable due to ZT_max being out of range
            break
        end
        
        % energy resolved quantities (for plotting ZT vs. b load line)
        b_resolved(kappa_l_index,:) = this_b;
        b_tot_resolved(kappa_l_index,:) = this_b_tot;
        zt_resolved(kappa_l_index,:) = this_zt;
        
        % obtain various quantities at max ZT (for plotting ZT vs. b)
        % due to the sensitivity of ZT vs. b, refinement is sometimes needed to obtain a smooth curve
        % interpolate to the maximum value via SPline instead of taking a simple maximum of existing data
        refined_e_grid = linspace(fermi_list(index_at_max_zt-1), fermi_list(index_at_max_zt+1),1000);
        
        refined_zt = interp1(fermi_list, this_zt, refined_e_grid,'spline');
        refined_b = interp1(fermi_list, this_b, refined_e_grid,'spline');
        refined_b_tot = interp1(fermi_list, this_b_tot, refined_e_grid,'spline');
        refined_s = interp1(fermi_list, seebeck_list, refined_e_grid, 'spline');
        refined_e_cond = interp1(fermi_list, e_cond_list, refined_e_grid, 'spline');
        refined_l = interp1(fermi_list, lorenz_list, refined_e_grid, 'spline');        
        
        [zt_at_max(kappa_l_index), refined_index_at_max_zt] = max(refined_zt);
        b_at_max(kappa_l_index) = refined_b(refined_index_at_max_zt);
        b_tot_at_max(kappa_l_index) = refined_b_tot(refined_index_at_max_zt);
        s_at_max(kappa_l_index) = refined_s(refined_index_at_max_zt);
        e_cond_at_max(kappa_l_index) = refined_e_cond(refined_index_at_max_zt);
        l_at_max(kappa_l_index) = refined_l(refined_index_at_max_zt);
        
    end
    
    % dynamically refine kappa_l values to cover entire b range
    b_within_range_list = (b_at_max > min(b_range) & b_at_max < max(b_range));
    b_within_range = sum(b_within_range_list);
    
    kappa_l_lower = kappa_l_list(find(b_within_range_list == 1, 1, 'first')-1); % new kappa_l lower limit
    kappa_l_upper = kappa_l_list(find(b_within_range_list == 1, 1, 'last')+1); % new kappa_l upper limit
end

% quantities at max ZT that can be calculated from others
kappa_l_at_max = T*kb^2.*e_cond_at_max./b_at_max; % (W/K-m) lattice thermal conductivity
kappa_e_at_max = l_at_max.*kb^2.*e_cond_at_max; % (W/K-m) electronic thermal conductivity

% Locate the measured data load lines in the ZT vs. b plot
zt_ll    = interp1(kappa_l_list, zt_resolved, target_kappa_l);
b_ll     = interp1(kappa_l_list, b_resolved, target_kappa_l);
b_tot_ll = interp1(kappa_l_list, b_tot_resolved, target_kappa_l);

%% Output plots
fig_num = 1;

%-----------------------------------
% Plot #1: ZT|_{max} vs. b_l
%-----------------------------------

figure(fig_num); fig_num = fig_num+1;

plot(b_at_max, zt_at_max,'-'); % ZT|_{max} vs. b_l
hold on
plot(b_ll, zt_ll,'--r');
plot(target_b, ZT,'o','MarkerFaceColor','red');

legend('ZT|_{max}', ['\kappa_{L}=' sprintf('%1.2f W/K-m',target_kappa_l)], 'Input')

xlabel('b_L')
ylabel('ZT|_{max}')
plot_setup(b_at_max, zt_at_max, ZT,'');

%-----------------------------------
% Plot #2: ZT|_{max} vs. b_total
%-----------------------------------

figure(fig_num); fig_num = fig_num+1;

plot(b_tot_at_max, zt_at_max,'-'); % ZT|_{max} vs. b_l
hold on
plot(b_tot_ll, zt_ll,'--r');
plot(target_b_tot, ZT,'o','MarkerFaceColor','red');

legend('ZT|_{max}', ['\kappa_{tot}=' sprintf('%1.2f',kappa)], 'Input')

xlabel('b_{total}')
ylabel('ZT|_{max}')
plot_setup(b_tot_at_max, zt_at_max, ZT,'');

%-----------------------------------
% Plot #3: Seebeck coefficient vs. b
%-----------------------------------

figure(fig_num); fig_num = fig_num+1;

plot(b_at_max, s_at_max.*1e6, '-')
hold on
plot(target_b, seebeck_coeff*1e6, 'o','MarkerFaceColor','red')

legend('@ZT|_{max}', 'Input')

xlabel('b_L')
ylabel('|Seebeck coefficient| (uV/k)')
plot_setup(b_at_max, s_at_max.*1e6, seebeck_coeff*1e6,'');

%-----------------------------------
% Plot #4: Electrical resistivity vs. b
%-----------------------------------

figure(fig_num); fig_num = fig_num+1;

semilogy(b_at_max, 1./e_cond_at_max.*1e2, '-')
hold on
plot(target_b, resistivity, 'o','MarkerFaceColor','red')

legend('@ZT|_{max}','Input')

xlabel('b_L')
ylabel('Resistivity (\Omega-cm)')
plot_setup(b_at_max, 1./e_cond_at_max.*1e2, resistivity,'exempt');

%-----------------------------------
% Plot #5: Thermal conductivity vs. b
%-----------------------------------

figure(fig_num); fig_num = fig_num+1;

semilogy(b_at_max, kappa_l_at_max, 'sb')
hold on
plot(b_at_max, kappa_e_at_max, '<r')
plot(b_at_max, kappa_l_at_max+kappa_e_at_max, 'og')

plot(target_b, target_kappa_l, 's','MarkerFaceColor','blue')
plot(target_b, e_kappa_at_seebeck, '<','MarkerFaceColor','red')
plot(target_b, kappa, 'o','MarkerFaceColor','green')

legend('kappa_l @ZT|_{max}','kappa_e @ZT|_{max}','kappa_{total} @ZT|_{max}','Input kappa_l','Input kappa_e','Input kappa_{total}')

xlabel('b_L')
ylabel('Thermal conductivity (W/K-m)')
plot_setup(b_at_max, kappa_l_at_max+kappa_e_at_max, kappa,'exempt');

%-----------------------------------
% Plot #6: Lorenz number vs. b
%-----------------------------------

figure(fig_num); fig_num = fig_num+1;

plot(b_at_max, l_at_max, '-')
hold on
plot(target_b, e_kappa_at_seebeck*resistivity*1e-2/(kb^2)/T,'o','MarkerFaceColor','red')

legend('L @ZT|_{max}','Input')

xlabel('b_L')
ylabel('Lorenz number (/(k_B/q)^2)')
plot_setup(b_at_max, l_at_max, e_kappa_at_seebeck*resistivity*1e-2/(kb^2)/T,'');

end % END of function zt_vs_b

function plot_setup(b_at_max, line_data, pt_data, plot_option)
% Special setup function for plots
font_size = 12;

% adjust xlim to cover entire b
xlim([0 max(b_at_max)])

% adjust ylim to cover input data
if ~strcmp(plot_option,'exempt')
    line_data_min_diff = abs(max([(pt_data - min(line_data)) pt_data]));
    line_data_max_diff = abs(pt_data - max(line_data));
    
    if line_data_min_diff/line_data_max_diff > 20
        % pt_data too close to the minimum value
        ylim([pt_data-line_data_max_diff*20 pt_data+line_data_max_diff])
    elseif line_data_max_diff/line_data_min_diff > 20
        % pt_data too close to the maximum value
        ylim([pt_data-line_data_min_diff pt_data+line_data_min_diff*20])
    end
end

% image size adjustment
pos = get(gcf,'Position');
set(gcf,'Units','inches', 'Position',[pos(1) pos(2) 4 3.125]);

ax = gca; ax.YColor = 'black';

% setup plot background and fonts
set(gcf, 'color', 'white');
set(gca, 'FontSize', font_size);
set(get(gca,'title'), 'FontSize',font_size);
set(get(gca,'xlabel'),'FontSize',font_size);
set(get(gca,'ylabel'),'FontSize',font_size);
end % END of function plot_setup

% END OF FILE
