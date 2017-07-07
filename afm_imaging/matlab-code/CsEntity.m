classdef (Abstract) CsEntity
    %CSENTITY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        xref_vec;
        yref_vec;
        index_vec;
    end
    
    methods (Abstract)
        asVector(obj)
    end
    
end

