function axial = generate_axial_37(ring_shift)
%GENERATE_AXIAL_37  Axial coords for 37-tile mosaic, with optional per-ring rotation.
%   axial = generate_axial_37()
%   axial = generate_axial_37([s1 s2 s3])   % CCW index shifts for rings 1..3
%   - s1 affects ring 1 (tiles 2..7, size 6)
%   - s2 affects ring 2 (tiles 8..19, size 12)
%   - s3 affects ring 3 (tiles 20..37, size 18)

    % ---- normalize ring_shift safely ----
    if nargin < 1 || isempty(ring_shift), ring_shift = [0 -1 -2]; end
    ring_shift = double(ring_shift(:).');      % row vector
    if numel(ring_shift) < 3, ring_shift(end+1:3) = 0; end
    if numel(ring_shift) > 3, ring_shift = ring_shift(1:3); end
    ring_shift = round(ring_shift);            % ensure integer shifts

    axial = zeros(37,2);

    % center
    axial(1,:) = [0 0];

    % axial directions for flat-top hex grid
    dirs = [ 1  0;   % E
             0  1;   % NE
            -1  1;   % NW
            -1  0;   % W
             0 -1;   % SW
             1 -1];  % SE

    % walk order used before: NW, W, SW, SE, E, NE
    order = [3 4 5 6 1 2];

    writePos = 2;
    ranges = {2:7, 8:19, 20:37};      % ring index ranges (sizes 6, 12, 18)

    for k = 1:3
        q = k; r = 0;                 % start at (k,0)
        tmp = zeros(6*k, 2);
        t=1;
        for d = order
            for s = 1:k
                tmp(t,:) = [q r];
                t = t + 1;
                q = q + dirs(d,1);
                r = r + dirs(d,2);
            end
        end
        % rotate ring k CCW by ring_shift(k) positions
        n = size(tmp,1);
        sh = mod(ring_shift(k), n);
        idx = mod((0:n-1) + sh, n) + 1;
        axial(ranges{k},:) = tmp(idx, :);
    end
end
