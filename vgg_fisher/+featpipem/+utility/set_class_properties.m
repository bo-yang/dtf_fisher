function rem_args = set_class_properties(obj, varargin)
% SET_CLASS_PROPERTIES  Parse list of parameter-value pairs
%   REM_ARGS = SET_CLASS_PROPERTIES(OBJ, ARGS) updates the object OBJ based on
%   the specified parameter-value pairs ARGS={PAR1, VAL1, ... PARN,
%   VALN}. The function produces an error if an unknown parameter name
%   is passed in.
%
%   Example::
%     The function can be used to parse a list of arguments
%     passed to a MATLAB functions:
%
%       function myFunction(x,y,z,varargin)
%       conf.parameterName = defaultValue;
%       set_class_properties(conf, varargin);

if ~isobject(obj), error('OBJ must be an object') ; end

remainingArgs = {} ;

args = varargin{:};

if mod(length(args),2) == 1
    error('Parameter-value pair expected (missing value?).') ;
end

for ai = 1:2:length(args)
    paramName = args{ai} ;
    value = args{ai+1} ;
    if isprop(obj,paramName)
        set(obj,paramName,value);
    else
        if nargout < 1
            error('Unknown parameter ''%s''.', paramName);
        else
            remainingArgs(end+1:end+2) = args(ai:ai+1);
        end
    end
end

rem_args = remainingArgs;
