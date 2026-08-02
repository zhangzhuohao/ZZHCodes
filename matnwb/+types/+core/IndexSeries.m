classdef IndexSeries < types.core.TimeSeries & types.untyped.GroupClass
% INDEXSERIES - Stores indices to image frames stored in an ImageSeries. The purpose of the ImageIndexSeries is to allow a static image stack to be stored somewhere, and the images in the stack to be referenced out-of-order. This can be for the display of individual images, or of movie segments (as a movie is simply a series of images). The data field stores the index of the frame in the referenced ImageSeries, and the timestamps array indicates when that image was displayed.
%
% Required Properties:
%  data, data_unit, indexed_timeseries


% REQUIRED PROPERTIES
properties
    indexed_timeseries; % REQUIRED ImageSeries
end

methods
    function obj = IndexSeries(varargin)
        % INDEXSERIES - Constructor for IndexSeries
        %
        % Syntax:
        %  indexSeries = types.core.INDEXSERIES() creates a IndexSeries object with unset property values.
        %
        %  indexSeries = types.core.INDEXSERIES(Name, Value) creates a IndexSeries object where one or more property values are specified using name-value pairs.
        %
        % Input Arguments (Name-Value Arguments):
        %  - comments (char) - Human-readable comments about the TimeSeries. This second descriptive field can be used to store additional information, or descriptive information if the primary description field is populated with a computer-readable string.
        %
        %  - control (uint8) - Numerical labels that apply to each time point in data for the purpose of querying and slicing data by these values. If present, the length of this array should be the same size as the first dimension of data.
        %
        %  - control_description (char) - Description of each control value. Must be present if control is present. If present, control_description[0] should describe time points where control == 0.
        %
        %  - data (int32) - Index of the frame in the referenced ImageSeries.
        %
        %  - data_continuity (char) - Optionally describe the continuity of the data. Can be "continuous", "instantaneous", or "step". For example, a voltage trace would be "continuous", because samples are recorded from a continuous process. An array of lick times would be "instantaneous", because the data represents distinct moments in time. Times of image presentations would be "step" because the picture remains the same until the next timepoint. This field is optional, but is useful in providing information about the underlying data. It may inform the way this data is interpreted, the way it is visualized, and what analysis methods are applicable.
        %
        %  - data_conversion (single) - Scalar to multiply each element in data to convert it to the specified 'unit'. If the data are stored in acquisition system units or other units that require a conversion to be interpretable, multiply the data by 'conversion' to convert the data to the specified 'unit'. e.g. if the data acquisition system stores values in this object as signed 16-bit integers (int16 range -32,768 to 32,767) that correspond to a 5V range (-2.5V to 2.5V), and the data acquisition system gain is 8000X, then the 'conversion' multiplier to get from raw data acquisition values to recorded volts is 2.5/32768/8000 = 9.5367e-9.
        %
        %  - data_resolution (single) - Smallest meaningful difference between values in data, stored in the specified by unit, e.g., the change in value of the least significant bit, or a larger number if signal noise is known to be present. If unknown, use -1.0.
        %
        %  - data_unit (char) - Base unit of measurement for working with the data. Actual stored values are not necessarily stored in these units. To access the data in these units, multiply 'data' by 'conversion'.
        %
        %  - description (char) - Description of the time series.
        %
        %  - indexed_timeseries (ImageSeries) - Link to ImageSeries object containing images that are indexed.
        %
        %  - starting_time (double) - Timestamp of the first sample in seconds. When timestamps are uniformly spaced, the timestamp of the first sample can be specified and all subsequent ones calculated from the sampling rate attribute.
        %
        %  - starting_time_rate (single) - Sampling rate, in Hz.
        %
        %  - timestamps (double) - Timestamps for samples stored in data, in seconds, relative to the common experiment master-clock stored in NWBFile.timestamps_reference_time.
        %
        % Output Arguments:
        %  - indexSeries (types.core.IndexSeries) - A IndexSeries object
        
        obj = obj@types.core.TimeSeries(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'indexed_timeseries',[]);
        misc.parseSkipInvalidName(p, varargin);
        obj.indexed_timeseries = p.Results.indexed_timeseries;
        if strcmp(class(obj), 'types.core.IndexSeries')
            cellStringArguments = convertContainedStringsToChars(varargin(1:2:end));
            types.util.checkUnset(obj, unique(cellStringArguments));
        end
    end
    %% SETTERS
    function set.indexed_timeseries(obj, val)
        obj.indexed_timeseries = obj.validate_indexed_timeseries(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'int32', val);
        types.util.validateShape('data', {[Inf]}, val)
    end
    function val = validate_indexed_timeseries(obj, val)
        val = types.util.validateSoftLink('indexed_timeseries', val, 'types.core.ImageSeries');
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.TimeSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        refs = obj.indexed_timeseries.export(fid, [fullpath '/indexed_timeseries'], refs);
    end
end

end