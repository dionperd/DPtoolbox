function axesScale=TS2DPscale(ax)

%This function returns the axis' scale of  a tstool achse object given its resolution and delta value

axesScale = resolution(ax); %if scale is 'linear' or 'arbitrary' we are done

%if it is logarithmic we have to figure out the delta value
if strcmpi(axesScale,'logarithmic');
    
    axesDelta = delta(ax);
    
    if (axesDelta==2)
        axesScale ='log2';
    elseif (axesDelta==10)
        axesScale ='log10';
    elseif (axesDelta==exp(1))
        axesScale ='ln';
    else
        error('axis'' is none of ''log2'', ''log10'', or ''ln'', since axis'' delta~={2,10,e}')
    end
end

