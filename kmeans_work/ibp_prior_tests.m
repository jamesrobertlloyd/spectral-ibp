colormap(jet)
possmat = [ones(10,1); zeros(10,1)];
% do shift moves and return results
possizes = [];
numv = 2.^[1:size(possmat,1)];
%% do naive search and evaluation
itr = 1;
while ~isempty(possmat)
    itr
    itr = itr + 1;
    new_possmat = [];
    numhash = containers.Map(uint64(0),true);
    for col=1:size(possmat,2)
        % perform shift moves
        selcol = possmat(:,col);
        chinds = find(diff(selcol) == -1);
        for ii =1:length(chinds)
            ind = chinds(ii);
            tmp = selcol;
            tmp(ind,1) = 0;
            tmp(ind+1,1) = 1;
            tmp_num = numv*tmp;
            if ~numhash.isKey(tmp_num)
                new_possmat(:, end+1) = tmp;
                numhash(tmp_num) = true;
            end
        end
    end
    possmat = new_possmat;
    imsc(possmat,'colormap',jet)
    input('enter to continue');
    possizes(end+1) = size(possmat,2);
end

%% do recursive search