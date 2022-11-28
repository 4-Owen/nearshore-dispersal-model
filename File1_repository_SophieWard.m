% ---------PARALLEL PTM, AARON OWEN AND SOPHIE WARD --------------
% input needed
% 1. slf output file from TELEMAC
% 2. Resaved velocities in .mat format
% 3. particle release sites - .mat
%------------------------------------------------------------------------
close all, clear all, clc

% Platform specific settings
if isunix
    disp('Running on a Linux platform')
    project_dir = '/XXX/';
    toolboxes_dir = strcat(project_dir, 'toolboxes/');
    results_dir = strcat(project_dir, 'results/');
    figures_dir = strcat(project_dir, 'figures/');
    file = strcat(project_dir, '../XXX.slf');
    myFolder = strcat(project_dir, 'XXX/');
    % Set via slurm script
    site_idx = str2num(getenv('site_idx'));
    n_tides = str2num(getenv('n_tides'));
    release_idx = str2num(getenv('release_idx'));
    particle_idx_start = str2num(getenv('particle_idx_start'));
    particle_idx_end = str2num(getenv('particle_idx_end'));
    nparticles = particle_idx_end - particle_idx_start;
    plotfig = 0;

elseif ispc % Local development machine paths
    disp('Running on a Windows platform')
    project_dir = 'XXX/';
    toolboxes_dir = strcat(project_dir, 'toolboxes\');
    results_dir = strcat(project_dir, 'results\');
    figures_dir = strcat(project_dir, 'figures\');
    file = strcat(project_dir, '\output.slf');
    % if on linux these set via slurm script
    site_idx = 1;
    tide_idx = 2;
    release_idx = 24;
    particle_idx_start = 1;
    particle_idx_end = 1;
    nparticles = particle_idx_end - particle_idx_start;
    plotfig = 1;
else
    disp('Platform not supported')
end

% Paths and Toolboxes
addpath(strcat(toolboxes_dir, 'telemac_tools'))
addpath(strcat(toolboxes_dir, 'latlonutm'))
savedsites = strcat(project_dir, 'SiteFile.mat');

%% User defined variables
hmin = .1; % Set minimumn water depth for bouncing back from shallows
tr = 3600/1; % Decide release frequency
tdev = 500; % Number timesteps for development and testing, <nsteps
tolerance = 1000; % Buffer around the mask
in_mask_check = true; % Handle scenarios where a particle's coordinates matches that of a node in the mask
grid_x = 100; % number of divisions of domain for the search grid
grid_y = 100; % number of divisions of domain for the search grid

m = telheadr(file); % Read telemac header file
DT = m.DT;
utmE = m.XYZ(:, 1); % X coordinates
utmN = m.XYZ(:, 2); % Y coordinates
ikle = m.IKLE;
daily = 24 * (3600 / DT); % number timesteps per daily file

load(savedsites); %load sites
Xpsite = xsite;
Ypsite = ysite;
xstart = Xpsite;
ystart = Ypsite;

ndays = 90;
nrel_end = ndays * 24 * (3600 / DT); % set number days running
nsteps = nrel_end; % define size of the
ndaystart = XXX;

%% Adaptive Region Of Interest

% sort domain into search boxes
x_min = min(utmE) - tolerance;
x_max = max(utmE) + tolerance;
lx = linspace(x_min, x_max, grid_x + 1);
y_min = min(utmN) - tolerance;
y_max = max(utmN) + tolerance;
ly = linspace(y_min, y_max, grid_y + 1);

% Calculate bounding rectangle
x_offset = (x_max - x_min) / grid_x;
y_offset = (y_max - y_min) / grid_y;

nsites = 1;
ntide = 1;
dtsave = 0;

% Cache
file_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

% To take full advantage of the file cache we need to run multiple tides in a single slurm job
for tide_idx = 1:n_tides

    %% RESULTS NAME
    ptm_filename = sprintf('site-%d_tide-%d_release-%d_particle_start-%d_particle_end-%d.mat', ...
    site_idx, tide_idx, release_idx, particle_idx_start, particle_idx_end)
    ptm_filepath = strcat(results_dir, ptm_filename)

    % Figure name
    figure_filename = sprintf('site-%d_tide-%d_release-%d_particle_start-%d_particle_end-%d.fig', ...
    site_idx, tide_idx, release_idx, particle_idx_start, particle_idx_end)
    figure_filepath = strcat(figures_dir, figure_filename)

    close all

    if plotfig == 1
        figure(1)
        hFig = figure(1);
        set(gcf, 'PaperPositionMode', 'auto')
        set(hFig, 'Position', [0 0 1440 800])
        hold on
        scatter(utmE, utmN, 'x')
    end

    save_idx = 1;

    % Individual data slices will be written to file and collated in post-processing
    xsave = NaN(nparticles, nsteps); % like this get 0s in the vector
    ysave = NaN(nparticles, nsteps);
    hsave = NaN(nparticles, nsteps);
    hstart = NaN(nparticles);
    lonsave = NaN(nparticles, nsteps);
    latsave = NaN(nparticles, nsteps);
    dttotal = zeros(nparticles); % Count distance travelled			 - post-processing?
    dtcum = zeros(nparticles, nsteps); % Count cumlative distance travelled - post-processing?
    ndie = NaN(nparticles);

    for p = particle_idx_start:particle_idx_end

        x1 = Xpsite(site_idx);
        y1 = Ypsite(site_idx);
        roi_size = 0.2; % Region of interest - adaptive

        nfnew = 0

        for nnt = release_idx:nrel_end + release_idx

            nf = ceil(nnt / daily); %to decide which file it's reading from

            if nf > nfnew % avoiding loading files each time
                clear U V H
                nfnew = nf;

                % Check if filename_1 is already loaded into our file cache
                filename_1 = [myFolder, 'file_' num2str(nf + ndaystart - 1) '.mat'];

                if isKey(file_cache, filename_1)
                    % Use U,V,H from file cache, faster to read from memory
                    UU(:, ((nf * daily) - daily + 1):(nf * daily)) = file_cache(filename_1).U;
                    VV(:, ((nf * daily) - daily + 1):(nf * daily)) = file_cache(filename_1).V;
                    HH(:, ((nf * daily) - daily + 1):(nf * daily)) = file_cache(filename_1).H;
                else
                    % Otherwise, load in the file and add a new map entry to our file cache for future reads
                    load(filename_1);
                    UU(:, ((nf * daily) - daily + 1):(nf * daily)) = U;
                    VV(:, ((nf * daily) - daily + 1):(nf * daily)) = V;
                    HH(:, ((nf * daily) - daily + 1):(nf * daily)) = H;
                    cache_struct = struct();
                    cache_struct.U = U;
                    cache_struct.V = V;
                    cache_struct.H = H;
                    file_cache(filename_1) = cache_struct;
                    clear cache_struct U V H;
                end

                % Check if filename_2 is already loaded into our file cache
                nff = nf + 1; % need to load second day too, because use nnt+1 below
                filename_2 = [myFolder, 'file_' num2str(nff + ndaystart - 1) '.mat'];

                if isKey(file_cache, filename_2)
                    % Use U,V,H from file cache, faster to read from memory
                    UU(:, ((nff * daily) - daily + 1):(nff * daily) - daily + 6) = file_cache(filename_2).U(:, 1:6);
                    VV(:, ((nff * daily) - daily + 1):(nff * daily) - daily + 6) = file_cache(filename_2).V(:, 1:6);
                    HH(:, ((nff * daily) - daily + 1):(nff * daily) - daily + 6) = file_cache(filename_2).H(:, 1:6);
                else
                    % Otherwise, load in the file and add a new map entry to our file cache for future reads
                    load(filename_2);
                    UU(:, ((nff * daily) - daily + 1):(nff * daily) - daily + 6) = U(:, 1:6);
                    VV(:, ((nff * daily) - daily + 1):(nff * daily) - daily + 6) = V(:, 1:6);
                    HH(:, ((nff * daily) - daily + 1):(nff * daily) - daily + 6) = H(:, 1:6);
                    % Add new file cache entry
                    cache_struct = struct();
                    cache_struct.U = U;
                    cache_struct.V = V;
                    cache_struct.H = H;
                    file_cache(filename_2) = cache_struct;
                    clear cache_struct U V H;
                end

                U = UU;
                V = VV;
                H = HH;

            else % don't load anything, already loaded
            end

            if nnt == release_idx
                xold = x1; yold = y1;
                xsave(save_idx, nnt) = x1;
                ysave(save_idx, nnt) = y1;
            else
                xold = xnew; yold = ynew;
            end

            [Uq, Vq, Hq] = unstruc_interp_telemac(xold, yold, utmE, utmN, U, V, H, ikle, nnt, roi_size, x_offset, y_offset, in_mask_check); %plots
            xnew = xold + (Uq * DT);
            ynew = yold + (Vq * DT);
            [Unew, Vnew, Hnew] = unstruc_interp_telemac(xnew, ynew, utmE, utmN, U, V, H, ikle, nnt + 1, roi_size, x_offset, y_offset, in_mask_check);

            if Hnew < hmin
                xnew = xold;
                ynew = yold;
                ndie(p) = nnt;
            end

            b_threshold = 1000; %threshold if within this close of the boundary

            if ((ynew > (min(utmN) + b_threshold) && ynew < (max(utmN) - b_threshold))) && ((xnew > (min(utmE) + b_threshold) && xnew < (max(utmE) - b_threshold)))
                xsave(save_idx, nnt + 1) = xnew;
                ysave(save_idx, nnt + 1) = ynew;
                hsave(save_idx, nnt + 1) = Hnew;
                dtsave = dtsave + sqrt(((Uq * DT).^2) + (Vq * DT).^2);
                dttotal(save_idx) = dtsave; % Total dist
                dtcum(save_idx, nnt + 1) = dttotal(save_idx) + sqrt(((Uq * DT).^2) + (Vq * DT).^2); % Cumulative
                loopbreak = 0;
            else
                loopbreak = 1;
                break
            end

        end

        save_idx = save_idx + 1;
        clear U V H UU VV HH
    end %p=particle

    xsave(xsave == 0) = NaN;
    ysave(ysave == 0) = NaN;

    save(ptm_filepath, 'utmE', 'utmN', 'xstart', 'ystart', 'DT', 'xsave', 'ysave', 'hsave', '-v7.3')

    clear xsave ysave hsave hstart lonsave latsave dttotal dtcum ndie xold yold Uq Vq Hq Unew Vnew Hnew
end
