function delete_this(mesh,varargin)

narginchk(1,7);
nargoutchk(0,2);

% default
order = 2;
central = false;
second_mesh_basis = [];
solver = get_solver;
OPTIONS = solver_options;
isAll = false;

if ~isempty(varargin)
    if length(varargin) >= 1
        % frequency
        if ~ischar(varargin{1}) && ~isstring(varargin{1}) && (numel(varargin{1})==1)
            % sanity check
            if varargin{1} > 2
                error('Maximum supported moment order is 2. Please see help for details on how to use this function.')
            else
                order = 2; 
            end
        else
            error('Bad 2nd argument value. A scalar expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 2
        % do we want centralised moments
        if ischar(varargin{2}) || isstring(varargin{2})
            if strcmp(varargin{2},'central')
                central = true;
            end
        end
    end
    if length(varargin) >= 3
        % second mesh basis
        if isstruct(varargin{3}) || isempty(varargin{3})
            second_mesh_basis = varargin{3};
        else
            error('Bad 3nd argument value. A mesh structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 4
        % solver
        if ischar(varargin{4}) || isstring(varargin{4})
            % user specified, sanity check
            solver = get_solver(varargin{4});
        elseif isstruct(varargin{4})
            OPTIONS = varargin{4};
        elseif isempty(varargin{4})
            solver = get_solver;
        else
            error('Bad 4th argument value. Solver name or solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 5
        % solver options
        if isstruct(varargin{5})
            OPTIONS = varargin{5};
        elseif isempty(varargin{5})
            OPTIONS = solver_options;
        else
            error('Bad 5th argument value. Solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) == 6
        % if scattering for CW data as well
        if ischar(varargin{6}) || isstring(varargin{6})
            if strcmp(varargin{6},'all')
                isAll = true;
            end
        else
            error('Bad 6th argument value. Text expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) > 6
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

disp(['order ', num2str(order)])
disp(['central ', num2str(central)])
disp(['second_mesh_basis ', second_mesh_basis])
disp(['solver ', solver])
disp('options '); disp(OPTIONS)
disp(['all ', num2str(isAll)])

end

