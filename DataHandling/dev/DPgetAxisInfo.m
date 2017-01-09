function [axesLabels, axesUnits, axesScales, axesValues] = DPgetAxisInfo(axesCell)

%This function unloads the labels, units, scales and values of a cell of
%tstool achse objects

axesLabels=cellfun(@(ax) name(ax),axesCell,'uniformoutput',false);
axesUnits=cellfun(@(ax) char(unit(ax)),axesCell,'uniformoutput',false);
axesScales=cellfun(@(ax)TS2DPscale(ax),axesCell,'uniformoutput',false);
axesValues=cellfun(@(ax) values(ax),axesCell,'uniformoutput',false);


