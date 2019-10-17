function yhat = lsqisotonic(x,y,w)

%Code by Yun Wang

n = numel(x);

% Sort points ascending in x, break ties with y.
[xyord,ord] = sortrows([x(:) y(:)]); iord(ord) = 1:n;
yhat = xyord(:,2);

block = 1:n;
if (nargin == 3) && ~isempty(w)
    w = w(:); w = w(ord); % reorder w as a column

    % Merge zero-weight points with preceeding pos-weighted point (or
    % with the following pos-weighted point if at start).
    posWgts = (w > 0);
    if any(~posWgts)
        idx = cumsum(posWgts); idx(idx == 0) = 1;
        w = w(posWgts);
        yhat = yhat(posWgts);
        block = idx(block);
    end

else
    w = ones(size(yhat),class(yhat));
end
while true
    % If all blocks are monotonic, then we're done.
    diffs = diff(yhat);
    if all(diffs >= 0), break; end

    % Otherwise, merge blocks of non-increasing fitted values, and set the
    % fitted value within each block equal to a constant, the weighted mean
    % of values in that block.
    idx = cumsum([1; (diffs>0)]);
    sumyhat = accumarray(idx,w.*yhat);
    w = accumarray(idx,w);
    yhat = sumyhat ./ w;
    block = idx(block);
end

% Broadcast merged blocks out to original points, and put back in the
% original order and shape.
yhat = yhat(block);
yhat = reshape(yhat(iord), size(y));


end