function C = combine_STPC_smooth(A,B);

% Get field names of A and B
fieldsA = fieldnames(A);
fieldsB = fieldnames(B);

% Find common fields
commonFields = intersect(fieldsA, fieldsB);

% Check that A and B have the same number of elements
nA = numel(A);
nB = numel(B);

% Initialize C as a copy of A (so C already has all fields from A)
C = A;

% Loop over each common field name
for k = 1:numel(commonFields)
    f = commonFields{k};   % current field name

    % Concatenate values from A and B for this field
    % Here we use horizontal concatenation [ ... , ... ]
    % You can change to vertical concatenation [ ... ; ... ]
    % depending on your data shape.
    C(nA+1).(f) = B.(f);
    % If you want vertical concatenation instead, use:
    % C(i).(f) = [A(i).(f); B(i).(f)];
end

end

