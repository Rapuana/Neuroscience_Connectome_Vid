function [VEm, empty, emptyVals] = grabfiles(FB, PredAnti, Pers)
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
emptyVals = {}; % Create a cell array

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
                emptyVals = {emptyVals, char(f(ii))}; 
            else
                fprintf('There was an error! The message was:\n%s',...
                    e.message)
                fprintf('The identifier was:\n%s', e.identifier)
                error('Break loop and handle')
            end % error if
        end % try
    end % ismember if statement
end % f loop

emptyVals = emptyVals(2:end);

if nargout > 3
    error('Too many outputs defined')
end

    

end % function