function sidevid(outfilename, fps, timerange, Patient, FB, normalise)
%SIDEVID VEm1 is Predic, VEm2 is Anti
%
% AUTHOR: Samuel Dudley (dudley.physics@gmail.com) 2019

[VEm1, eA] = grabfiles(FB, 'Predic', Patient);
[VEm2, eP] = grabfiles(FB, 'Anti', Patient);

[m, nnodes] = size(VEm1);

if ~islogical(normalise)
    error('normalise must be true/false logical')
end

if m == size(VEm2,1) && nnodes == size(VEm2,2)
    % good
else
    error('VEm1 and VEm2 must be same size')
end

e = unique([eA, eP]);

for ii = 1:length(e)
    VEm1(:,e(ii)) = zeros(m,1);
    VEm2(:,e(ii)) = zeros(m,1);
end

if normalise
    VEm1STD = std(VEm1, 0, 1);
    VEm2STD = std(VEm2, 0, 1);
    VEm1mean = mean(VEm1, 1);
    VEm2mean = mean(VEm2, 1);
    
    for ii = 1:36
        VEm1(:,ii) = (VEm1(:,ii) - VEm1mean(ii) * ones(1898, 1)) ./ VEm1STD(ii);
        VEm2(:,ii) = (VEm2(:,ii) - VEm2mean(ii) * ones(1898, 1)) ./ VEm2STD(ii);
    end
end

load('MEGData.mat', 'time', 'Hz', 'MNI')

loc = MNI;

sphereradius = 15;

writerObj = VideoWriter(outfilename, 'MPEG-4');
writerObj.FrameRate = fps;
writerObj.Quality = 100; % best quality
open(writerObj)

outwidth = 1024;
outheight = 576;
hf = figure(888);
figpos = [50 50 outwidth outheight];

set(hf, 'position', figpos, 'color', [1 1 1], 'paperPositionMode', 'auto');


h1 = subplot(1,2,1);
set(h1, 'visible', 'off', ...
    'units', 'normalized', ...
    'position', [0 0 0.5 1], ...
    'color', 'w');

% ax1 = axes(h1, 'visible', 'off', ...
%     'units', 'normalized', ...
%     'position', [0 0 1 1], ...
%     'color', 'w');
colormap('jet')

crange = quantile([VEm1(:); VEm2(:)], [1e-3, 1-1e-3]);
caxis(crange);
cbar = colorbar;
text(0.5, 0.05, ...
    '\textbf{Predictable}', ...
    'units', 'normalized', ...
    'FontSize', 24, ...
    'HorizontalAlignment', 'center', ...
    'Interpreter', 'LaTeX');


set(h1, 'position', [0 0 0.5 1])
set(cbar, 'position', [0.47, 0.05, 0.02, 0.8],...
    'FontSize', 14, ...
    'TickLabelInterpreter', 'LaTeX')

h2 = subplot(1,2,2);
set(h2, 'visible', 'off', ...
    'units', 'normalized', ...
    'position', [0.5 0 0.5 1], ...
    'color', 'w');
text(0.5, 0.05, ...
    '\textbf{Anti}', ...
    'units', 'normalized', ...
    'FontSize', 24, ...
    'HorizontalAlignment', 'center', ...
    'Interpreter', 'LaTeX');

text(0, 0.97, ...
    {['\textbf{Patient} ', num2str(Patient)]; FB}, ...
    'units', 'normalized', ...
    'FontSize', 24, ...
    'HorizontalAlignment', 'center', ...
    'Interpreter', 'LaTeX');


% VIEW AND CAMROLL




timerange = timerange * Hz / 1e3; % convert to loop values

timerange(1) = floor(timerange(1)); timerange(2) = ceil(timerange(2));

if timerange(1) < 1
    timerange(1) = 1;
end

if timerange(2) > time(end)
    timerange(2) = time(end);
end

for jj = timerange(1):timerange(2)
    if ~mod(jj, round(fps)*60)
        fprintf(1, '%d ', jj);
        if ~mod(jj, 20*round(fps)*60)
            fprintf(1, '\n')
        end
    end
    % Set up the spheres
    if jj == timerange(1)
        ylabel(cbar, 'Units [unknown]', ...
            'Interpreter', 'LaTeX', 'FontSize', 20)
        sp1 = add_sphere_size_internal(h1, loc(:,1), ...
                                           loc(:,2), ...
                                           loc(:,3), ...
                                           sphereradius*ones(nnodes,1), ...
                                           VEm1(1,:));
        sp2 = add_sphere_size_internal(h2, loc(:,1), ...
                                           loc(:,2), ...
                                           loc(:,3), ...
                                           sphereradius*ones(nnodes,1), ...
                                           VEm2(1,:));
        cdata1 = get(sp1, 'cdata');
        cdata2 = get(sp2, 'cdata');
        
        th = text(0.83, 0.95, ...
                sprintf('t = %.1f ms\n %s', time(jj), 'stimulus'), ...
                'units', 'normalized', ...
                'FontSize', 18, ...
                'HorizontalAlignment', 'center', ...
                'Interpreter', 'LaTeX');
    else
        for k = 1:nnodes
            cdata1(:, (k-1) * 22 + (1:21)) = VEm1(jj, k);
            cdata2(:, (k-1) * 22 + (1:21)) = VEm2(jj, k);
        end
        set(sp1, 'cdata', cdata1)
        set(sp2, 'cdata', cdata2)
        
        if time(jj) < 200
            stim = 'Pre-Stimulus';
        elseif time(jj) < 700
            stim = 'Stimulus 1 (200 ms)';
        elseif time(jj) < 1200
            stim = 'Stimulus 2 (700 ms)';
        elseif time(jj) < 1700
            stim = 'Stimulus 3 (1200 ms)';
        elseif time(jj) < 2200
            stim = 'Stimulus 4 (1700 ms)';
        else
            stim = 'Stimulus 5 (2200 ms)';
        end
            
        set(th, 'string', sprintf('t = %.1f ms\n%s', time(jj), stim))
    end
    
    writeVideo(writerObj, getframe(hf));
    
end

close(writerObj);

end % End function

function hout=add_sphere_size_internal(ax, x, y, z, s, c)
% input:    x   Centrepoint of the sphere in the x direction
%           y   Centrepoint of the sphere in the y direction
%           z   Centrepoint of the sphere in the z direction
%           s   Size for each sphere
%           c   Value by which to color the sphere
%           ax  Axis upon which we place the spheres.
%
% output: (optional) handles to the surfs
% 
%
% Anton Lord, UQ 2010
% James Roberts, QIMR Berghofer, 2014-2017

hold(ax,'on');

n = 20; % number of faces each sphere has
c = double(c);
[x0,y0,z0] = sphere(n);

nspheres=length(x);
Xall=nan(n+1,nspheres*(n+2)); % each sphere is an (n+1)-by-(n+1) matrix
Yall=Xall;
Zall=Xall;
Call=Xall;

for j = 1:length(x)
    if size(c,1) == 1
        intensity = zeros(n+1)+c(j); 
    end
    Xall(:,(j-1)*(n+2)+(1:(n+1)))=x0*s(j)+x(j);
    Yall(:,(j-1)*(n+2)+(1:(n+1)))=y0*s(j)+y(j);
    Zall(:,(j-1)*(n+2)+(1:(n+1)))=z0*s(j)+z(j);
    Call(:,(j-1)*(n+2)+(1:(n+1)))=intensity;
end
h=surf(ax,Xall,Yall,Zall,Call,'EdgeColor','none');
%light('Position',[15 -20 10])
daspect(ax,[1 1 1])
if nargout>0
    hout=h;
end
end

function [VEm, empty] = grabfiles(FB, PredAnti, Pers)
%GRABFILES [VEm, empty, emptyVals] = grabfiles(FB, PredAnti, Pers) is a 
% function used to obtain the values from the file locations. Here we look 

% FB        -> part1 = 'Body'; % Face or Body
% PredAnti  -> part2 = 'Predic'; % Anti or Predic 
% Pers      -> no = 6; % person

% Check inputs
if ischar(FB) && contains(FB, {'Face', 'Body'})
    % good
elseif isstring(FB) 
    FB = char(FB);
    if ~contains(FB, {'Face', 'Body'})
        error("Incorrect input of FB, it must contain 'Face' or 'Body'")
    end
    fprintf('Succesfully converted FB to character array from string.\n')
elseif ~contains(PredAnti, {'Face', 'Body'})
    fprintf("Expecting Character input for FB of either 'Face' or 'Body'\n")
    fprintf("Instead input is FB = '%s'\n", FB)
    error('Error')
else
    fprintf('Invalid class type for FB, expecting char array, class = %s', ...
        class(FB))
    error('Invalid input for FB')
end

if ischar(PredAnti) && contains(PredAnti, {'Anti', 'Predic'})
    % good
elseif isstring(PredAnti) 
    PredAnti = char(PredAnti);
    if ~contains(PredAnti, {'Anti', 'Predic'})
        error("Incorrect input of PredAnti, it must contain 'Predic' or 'Anti'")
    end
    fprintf('Succesfully converted PredAnti to character array from string.\n')
elseif ~contains(PredAnti, {'Predic', 'Anti'})
    fprintf("Expecting Character input for PredAnti of either 'Predic' or 'Anti'\n")
    fprintf("Instead input is PredAnti = '%s'\n", PredAnti)
    error('Error')
else
    fprintf('Invalid class type for PredAnti, expecting char array, class = %s', ...
        class(PredAnti))
    error('Invalid input for PredAnti')
end

if isnumeric(Pers)
    if Pers ~= floor(Pers)
        error('Expecting integer value for Pers')
    end
else
    error('Expecting an numeric integer value for Pers')
end

% Set up file locations
cd '/Users/Dudley/Documents/University/MXB371/Data';
Dir = cd;

files = dir('RawData');
f = {files.name};


% Now find hidden files beginning with or containing '.'
nonfiles = find(contains(f, '.'));

% Define length
nodes = length(f) - length(nonfiles);



% Now I need to get all the data. 

% Preset a matrix that will contain all VEm values. 
VEm = zeros(1898, nodes);

% set folder location
folDir = [Dir, filesep, 'RawData'];

% We want to save the files that are empty.
empty = [];

% iter
iter = 0;
for ii = 1:length(f)
    if ~ismember(ii,nonfiles)
        iter = iter + 1;
        % Firstly I need to read the data
        % Check file exists.
        filecheck = exist([folDir, filesep, char(f(ii))], 'file');
        if filecheck == 0
            warning('File does not exist')
        elseif filecheck == 7
            % good, ignore
        else
            error('Non desired file type already exists')
        end
        
        % Set the file location
        fileloc = [folDir, filesep, char(f(ii))];
        
        % Read in the files at this location
        b = dir(fileloc);
        b = {b.name};
        
        % Now we check how many files exist in here. 
        try
            load([fileloc, filesep, FB, PredAnti, '5_', num2str(Pers)], ...
                'orvem')
            VEm(:, iter) = orvem';
        catch e
            if strcmp(e.identifier, 'MATLAB:load:couldNotReadFile')
                fprintf('No such file exists for %s\n', char(f(ii)))
                fprintf('Count: %g\n f: %g\n', iter, ii)
                empty = [empty, iter]; % Just the way it gotta be lord.
            else
                fprintf('There was an error! The message was:\n%s',...
                    e.message)
                fprintf('The identifier was:\n%s', e.identifier)
                error('Break loop and handle')
            end % error if
        end % try
    end % ismember if statement
end % f loop

end % function grabfiles
