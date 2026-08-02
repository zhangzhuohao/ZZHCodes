classdef Images < types.core.NWBDataInterface & types.untyped.GroupClass
% IMAGES - A collection of images.
%
% Required Properties:
%  description, image


% REQUIRED PROPERTIES
properties
    description; % REQUIRED (char) Description of this collection of images.
    image; % REQUIRED (Image) Images stored in this collection.
end

methods
    function obj = Images(varargin)
        % IMAGES - Constructor for Images
        %
        % Syntax:
        %  images = types.core.IMAGES() creates a Images object with unset property values.
        %
        %  images = types.core.IMAGES(Name, Value) creates a Images object where one or more property values are specified using name-value pairs.
        %
        % Input Arguments (Name-Value Arguments):
        %  - description (char) - Description of this collection of images.
        %
        %  - image (Image) - Images stored in this collection.
        %
        % Output Arguments:
        %  - images (types.core.Images) - A Images object
        
        obj = obj@types.core.NWBDataInterface(varargin{:});
        [obj.image, ivarargin] = types.util.parseConstrained(obj,'image', 'types.core.Image', varargin{:});
        varargin(ivarargin) = [];
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        misc.parseSkipInvalidName(p, varargin);
        obj.description = p.Results.description;
        if strcmp(class(obj), 'types.core.Images')
            cellStringArguments = convertContainedStringsToChars(varargin(1:2:end));
            types.util.checkUnset(obj, unique(cellStringArguments));
        end
    end
    %% SETTERS
    function set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function set.image(obj, val)
        obj.image = obj.validate_image(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
        types.util.validateShape('description', {[1]}, val)
    end
    function val = validate_image(obj, val)
        constrained = { 'types.core.Image' };
        types.util.checkSet('image', struct(), constrained, val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBDataInterface(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        io.writeAttribute(fid, [fullpath '/description'], obj.description);
        refs = obj.image.export(fid, fullpath, refs);
    end
end

end