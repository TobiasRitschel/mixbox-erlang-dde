function [A, B, C] = createSystemMatrices(c, a, M, varargin)

% Return sparse matrices?
SparseMatrices = false;
if(nargin > 3 && strcmpi(varargin{1}, 'sparse')), SparseMatrices = true; end

% Number of delayed quantities
nz = numel(M);

switch SparseMatrices
    case true
        % Rows, columns, and values
        Ar = nan(sum(2*M + 1), 1);
        Ac = nan(sum(2*M + 1), 1);
        Av = nan(sum(2*M + 1), 1);

        Br = nan(nz, 1);
        Bc = nan(nz, 1);
        Bv = nan(nz, 1);

        Cr = nan(sum(M+1), 1);
        Cc = nan(sum(M+1), 1);
        Cv = nan(sum(M+1), 1);

        % Matrix dimensions
        Arows = sum(M) + nz;
        Acols = Arows;

        Brows = Arows;
        Bcols = nz;

        Crows = nz;
        Ccols = Acols;

        % Offsets
        offsetA = 0;
        offsetB = 0;
        offsetC = 0;
        offsetM = 0;

        for i = 1:nz
            % Parameters for i'th kernel
            Mi = M(i);
            ai = a(i);
            ci = c(i, 1:Mi+1);

            % Diagonal of the A matrix
            idx = offsetA + (1:Mi+1);
            Ar(idx) = offsetM + (1:Mi+1);
            Ac(idx) = offsetM + (1:Mi+1);
            Av(idx) = -ai;

            offsetA = offsetA + Mi+1;

            % Lower bidiagonal of the A matrix
            idx = offsetA + (1:Mi);
            Ar(idx) = offsetM + (2:Mi+1);
            Ac(idx) = offsetM + (1:Mi  );
            Av(idx) = ai;

            offsetA = offsetA + Mi;

            % The B matrix
            idx = offsetB + 1;
            Br(idx) = offsetM + 1;
            Bc(idx) = i;
            Bv(idx) = ai;

            offsetB = offsetB + 1;

            % The C matrix
            idx = offsetC + (1:Mi+1);
            Cr(idx) = repelem(i, Mi+1);
            Cc(idx) = offsetM + (1:Mi+1);
            Cv(idx) = ci;

            offsetC = offsetC + Mi+1;

            % Increment offset
            offsetM = offsetM + Mi+1;
        end

        % Create sparse matrix
        A = sparse(Ar, Ac, Av, Arows, Acols);
        B = sparse(Br, Bc, Bv, Brows, Bcols);
        C = sparse(Cr, Cc, Cv, Crows, Ccols);

    case false
        % Allocate memory
        Ai = cell(1, nz);
        bi = cell(1, nz);
        ci = cell(1, nz);

        for i = 1:nz
            % Unit vector
            e1 = eye(M(i)+1, 1);

            % Identity matrices
            IM = eye(M(i)+1, M(i)+1);

            % Lower shift matrix
            LM = diag(ones(M(i), 1), -1);

            % System matrices
            Ai{i} = a(i)*(LM - IM);
            bi{i} = a(i)* e1;
            ci{i} = c(i, 1:M(i)+1);
        end

        % Create block-diagonal matrices
        A = blkdiag(Ai{:});
        B = blkdiag(bi{:});
        C = blkdiag(ci{:});
end