classdef frdExt < frd
    %FRFEXT Summary of this class goes here
    %   Detailed explanation goes here
    % 
    %    frdExt(respData, freqs, Ts).
    % Calling sequence is exactly the same as frd
    % 
    % You can also include coherence data:
    %   frdExt(respData, freqs, CoherenceData, Ts)
    % Methods:
    %   frfBode:  overloads standard frfBode function. Always plots into
    %   units of Hz.
    %   [hMag, hPhase] = frfBode(obj, F1, plotstyle)
    %
    % See Also: frd, frfBode
    properties
        Coherence;
                
    end
    
    methods
        function obj = frdExt(varargin)
             if length(varargin{3}) == length(varargin{1})
                % The only reason the third argument should be "long", ie
                % same length as the response data is if the user included
                % coherence infor
                Coherence = varargin{3};
                varargin(3) = [];
            else
                Coherence = [];
             end
            
            [Gz, w_s, ind_s] = monotonicFRF(varargin{1}, varargin{2});
            obj@frd(Gz, w_s, varargin{3:end});
            if ~isempty(Coherence); Coherence(ind_s) = []; end;
            obj.Coherence = Coherence;
        end
        
        function PlotHandles = frfBode(obj, F1, plotstyle)
            if ~exist('F1', 'var') | isempty('F1')
                F1 = figure(100)
            end
            if ~exist('plotstyle', 'var') | isempty('plotstyle')
                plotstyle = 'b';
            end
            
%             [hMag, hPhase] = frfBode(squeeze(obj.ResponseData), obj.Frequency, F1,...
%                              plotstyle, obj.FrequencyUnit)
            tmpSys = chgFreqUnit(obj, 'Hz');
            if isempty(obj.Coherence)
                cohFlag = 0;
            else
                cohFlag = 1;
            end
            Subs = [211;212];
            PlotHandles = frfBode_internal(tmpSys,cohFlag, F1,Subs, plotstyle);

            
        end
        function H = interpFRF(obj, w_s)
            sz = size(obj.ResponseData);
            no = sz(2);
            ni = sz(1);
            for j=1:ni
                for i=1:no
                    H_frf(j,i,:) = interp1(obj.Frequency, squeeze(obj.ResponseData(j,i,:)), w_s, 'spline');
                end
            end
            H = frd(H_frf, w_s, obj.Ts);
            
            
        end
        
    end
    
end


function [PlotHandles] = frfBode_internal(SYS, cohFlag, F1,subs, plotStyle)

Gz_frf = squeeze(SYS.ResponseData(1,1,:));
w_s    = SYS.Frequency;

if ~exist('F1', 'var') || isempty(F1)
    F1 = figure;
end

if ~exist('plotStyle', 'var') || isempty(plotStyle)
    plotStyle = 'b';
end




magsDB = 20*log10(abs(Gz_frf));
phase  = unwrap(angle(Gz_frf))*180/pi;

figure(F1)
subplot(subs(1,:))
    if ~cohFlag
        [h_mag] = semilogx(w_s, magsDB, plotStyle, 'LineWidth', 1.5);
        PlotHandles.h_mag = h_mag;
        grid on; hold on;
        xlim([w_s(1), w_s(end)]);
    else
        
        [AX, h_mag, h_coh] = plotyy(w_s, magsDB, w_s, SYS.Coherence, 'semilogx', 'semilogx');

        grid on; hold on;
        AX(1).XLim = [w_s(1), w_s(end)];
        AX(2).XLim = [w_s(1), w_s(end)];
%         AX(2).YLim = [0,2];
        AX(1).YLim = [AX(1).YLim(1), max(magsDB)*1.1];
        [style, color] = splitLineStyle(plotStyle);
        if ~isempty(style)
        h_mag.LineStyle = style;
        end
        if ~isempty(color)
           h_mag.Color = color;
           AX(1).YColor = color;
        end

        ylabel(AX(1), 'Mag [dB]');
        ylabel(AX(2), 'Coherence');
        PlotHandles.h_mag = h_mag;
        PlotHandles.h_coh = h_coh;
        PlotHandles.AX    = AX;
        
        
    end
subplot(subs(2,:))
    h_phase = semilogx(w_s, phase, plotStyle, 'LineWidth', 1.5);
    grid on; hold on;
    xlim([w_s(1), w_s(end)]);
% if exist('freqUnit', 'var')
    xlabel(SYS.FrequencyUnit)
% end
    PlotHandles.h_phase = h_phase;
    
    ylabel('Phase [deg]')

end



function [style, color] = splitLineStyle(st)
    inds = isletter(st);
    if inds(end) == 1
        % contains color
        color = st(end);
        style = st(1:end-1);
    else
        % no color specified
       style = st;
       color = [];
    end
    


end




