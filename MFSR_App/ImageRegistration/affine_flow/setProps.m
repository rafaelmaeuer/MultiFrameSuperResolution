function obj = setProps(obj, varargin)
%setProps sets object properties
%   OBJ = setProps(OBJ, NAME1, VAL1, NAME2, VAL2, ...) sets property NAME1
%   of object OBJ to VAL1, property NAME2 to VAL2 and so on. The NAMEs must
%   be strings giving the names of properties with public set access in the
%   object.

% Copyright David Young 2010

% Line below isn't needed if this function is used as is. However, if it is
% copied into the methods section of a class, and if updating of private
% properties using this mechanism is not to be allowed, then the line below
% and the 3 lines further down must be reinstated
% pnames = propNames(obj, true);

len = length(varargin);
for inp = 1:2:len
    if inp == len
        error('setProps:noval', 'Missing value for final property');
    end
    prop = varargin{inp};
%     See note above: reinstate lines below if necessary
%     if ~ismember(prop, pnames)
%         error('setProps:badname', 'Unrecognised property name');
%     end
    obj.(prop) = varargin{inp+1};
end

end
