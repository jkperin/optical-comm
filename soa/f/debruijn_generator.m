function sequence = debruijn_generator(num_chars,subsequence_length)
%DEBRUIJN_GENERATOR compute example de Bruijn sequence
%   sequence = debruijn_generator(num_chars,subsequence_length)
%   where num_chars is an integer specifying the number of arbitrary
%   characters and subsequence_length is an integer specifying 
%   subsequence length. 
% 
%   This function will compute an example of an interleaved pseudo-random
%   de Bruijn sequence. The sequence is, by nature, circular, so the final
%   members of pairs (or triplets, quadruplets, etc...) are found wrapped
%   around to the beginning of the sequence.
%   Example:
%   >> sequence = debruijn_generator(3,2)'
%   sequence = 
%       [1 2 1 3 3 2 2 3 1]
%
%   A sequence of N characters with a subsequence length of L would be
%   ---------------------------------------------
%   N^L characters in length,and would contain:
%   N^(L-1) examples of each character,
%   N^(L-2) examples of each pair of characters,
%   N^(L-3) examples of each possible triplet of characters,
%   and so on.
%   ---------------------------------------------
%
%   Computation times will be *considerable* and unpredictable for large
%   numbers of characters (>10) and subsequence lengths (>4). Sometimes the
%   code will be unable to compute a solution and will restart. This
%   restart is initiated when the total computation time exceeds time taken 
%   before the first backtrack times a multiplier of 4 (arrived at
%   empirically).
%
%   It is advisable to reseed the random number generator in Matlab at the
%   start of each session.
% 
%   Author: W. Owen Brimijoin
%   MRC Institute of Hearing Research (Scottish Section)
%   owen(at)ihr.gla.ac.uk
%
%   If you find this function useful and use it for experimental design,
%   please consider citing Brimijoin and O'Neill, 2010 - (Hearing
%   Research).
%
%   Created: 30th November, 2008
%   Revised: 25th May, 2010 - Help text edited
%   Revised: 14th July, 2010 - Changed method for permuting all possible
%   combinations to that created by Jos van der Geest (combn.m available on
%   File Exchange). 

if nargin~=2
    display('Error: num_chars and subsequence_length are undefined')
        sequence = [];
    return
end

if floor(num_chars)~=num_chars | num_chars<1,
    display('Error: Number of characters must be positive integer')
    sequence = [];
    return
end

if floor(subsequence_length)~=subsequence_length | subsequence_length<1,
    display('Error: Subsequence length must be positive integer')
    sequence = [];
    return
end

if subsequence_length == 1,
    sequence = randperm(num_chars)';
    return
else
    % method for generating all possible combinations is from combn 
    % (c) Jos van der Geest:
    % create a list of all possible combinations of num_chars elements
    [combination_list{subsequence_length:-1:1}] = ndgrid(1:num_chars);
    % concatenate into one matrix, reshape into 2D and flip columns
    combination_list = reshape(cat(subsequence_length+1,combination_list{:}),[],subsequence_length);
end

sequence = [];
%this loop allows multiple attempts to compute the sequence. If an attempt
%has gone on too long (currently, first backtrack time * 4), the program
%will start again.
%==========================================================================
while isempty(sequence), 
    sequence = attempt_sequence(num_chars,subsequence_length,combination_list);
end
%==========================================================================


function sequence = attempt_sequence(num_chars,subsequence_length,combination_list)

sequence_length = num_chars^subsequence_length;
tic;
% create a progress display bar:
% progdisp = figure(...
%     'Units',            'normalized',...
%     'Position',         [0.5 0.5 0.25 0.02],...
%     'NumberTitle',      'off',...
%     'Resize',           'off',...
%     'MenuBar',          'none',...
%     'BackingStore',     'off' );
% progaxes = axes(...
%     'Position',         [0.02 0.15 0.96 0.70],...
%     'XLim',             [0 1],...
%     'YLim',             [0 1],...
%     'Box',              'on',...
%     'ytick',            [],...
%     'xtick',            [] );
% progpatch = patch(...
%     'XData',            [0 0 0 0],...
%     'YData',            [0 0 1 1],...
%     'EraseMode',        'none' );

sequence = []; %clear previous sequences
history = []; %clear previous sequence search history
backtrack_indices = []; %clear search route
%pick a starting group of characters at random:
history = randperm(size(combination_list,1));
backtrack_indices{1} = {history};
history = history(1); %add this choice to the search history

%use this starting group as the first element in the sequence:
sequence = combination_list(history,:)';
%current_spot tracks where the function is in the sequence generation:
current_spot = subsequence_length-(subsequence_length-2);
% progress = (length(sequence)/(sequence_length+2)); %update progress bar
% set(progpatch,'XData',[0 progress progress 0])
% set(progpatch,'FaceColor',[.55 .6 .8]);
% set(progdisp,'Name','Sequencing Permutations...')
% drawnow

%continue until sequence is complete...
backtrack_warning = 0;
while length(sequence)<sequence_length+(subsequence_length-1),
    match = sequence(end-(subsequence_length-2):end); %find the last bit of the sequence so far
    %find subsequences whose 1st members match last members of sequence
    working_combination_list = combination_list; %define an alternate character array
    working_combination_list(history,:) = 0; %zero those entries already chosen
    index = sum(abs(working_combination_list(:,1:subsequence_length-1)...
    - repmat(match',size(working_combination_list,1),1)),2); 
    index = find(index==0);
    backtrack_indices = [backtrack_indices;{index}];
    if sum(index)>0,
        %choose at random from matching subsequences:
        random_pick = randperm(length(index));
        random_pick = random_pick(1);
        index = index(random_pick); 
        %remove this choice from options:
        temp_backtrack = (backtrack_indices{end});
        temp_backtrack(temp_backtrack==index) = [];
        backtrack_indices{end} = {temp_backtrack};
        sequence = [sequence;combination_list(index,end)]; %add new subsequence to sequence
%         progress = (length(sequence)/(sequence_length+2)); %update progress bar
%         if ~ishandle(progdisp)
%             sequence = 'aborted';
%             return
%         end
%         set(progpatch,'XData',[0 progress progress 0])
%         set(progpatch,'FaceColor',[.55 .6 .8]);
%         if backtrack_warning == 0,
%             set(progdisp,'Name','Sequencing Permutations...')
%         else
%             set(progdisp,'Name','Backtracking and Choosing Alternates...')
%         end
%         drawnow
        history = [history;index];
    else
        if backtrack_warning == 0,
            attempt_duration = toc*4; %define the cutoff time for restart
        end
        backtrack_warning = 1;
        while length(cell2mat(backtrack_indices{end}))<1,
            if isempty(sequence), % if the code has backtracked to the beginning, 
                sequence = []; %clear sequence
                history = []; %clear previous sequence search history
                backtrack_indices = []; %clear search route
                %pick a starting group of characters at random:
                history = randperm(size(combination_list,1));
                backtrack_indices{1} = {history};
                history = history(1); %add this choice to the search history
                %use this starting group as the first element in the sequence:
                sequence = combination_list(history,:)';
                %current_spot tracks where the function is in the sequence generation:
                current_spot = subsequence_length-(subsequence_length-2);
            else % or if not,
                history(end) = []; %remove the last history entry
                backtrack_indices(end) = []; % remove last backtrack index
                sequence(end) = []; %remove last sequence character
            end
        end
        new_indices = cell2mat(backtrack_indices{end}); % look at the last backtrack index
        random_pick = randperm(length(new_indices)); % choose one at random
        index = new_indices(random_pick(1)); 
        temp_backtrack = cell2mat((backtrack_indices{end}));
        temp_backtrack(temp_backtrack==index) = []; %remove this choice from options
        backtrack_indices{end} = {temp_backtrack}; %add this choice to backtrack index
        sequence = [sequence;combination_list(index,end)]; %add new subsequence to sequence
%         progress = (length(sequence)/(sequence_length+2)); %update progress bar
%         if ~ishandle(progdisp)
%             sequence = 'aborted';
%             return
%         end
%         set(progpatch,'XData',[0 progress progress 0])
%         set(progpatch,'FaceColor',[.55 .6 .8]);
%         if backtrack_warning == 0,
%             set(progdisp,'Name','Sequencing Permutations...')
%         else
%             set(progdisp,'Name','Backtracking and Choosing Alternates...')
%         end
%         drawnow
        history = [history;index];
        if toc>attempt_duration %restart if the attempt has gone on too long.
            sequence = []; %empty out the sequence.
%             close(progdisp) %close progress bar
%             return
        end
    end

end

% close(progdisp) %close progress bar
sequence = sequence(1:sequence_length); %remove circular portion of sequence

%the end.