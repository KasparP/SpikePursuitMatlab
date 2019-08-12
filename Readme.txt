This is the version of the Spike Pursuit Code used to analyze data in March 2019, for revisions of Abdelfattah et al 2019 (Science)

This snapshot of the code is functional but is not intended as a practical tool for wide distribution; Spike Pursuit has since been translated to Python and integrated into CaIMAN by Andrea Giovannucci and lab - a development version may be available from Andrea.

Example data to accompany this code will soon be available via FigShare, and will be linked in this readme.

For a basic description of how Spike Pursuit works, see the supplement of Abdelfattah et al. 2019

The core function is spikePursuit.m, which calls denoiseSpikes.m in an inner loop.

spikePursuit.m loads data that has already been aligned and cut into manageable blocks (see an example 'data block' in the exampledata folder). The format for these blocks is very simple.

The alignment and block-cutting happens separately in the script forAmrita_loadData.

Spike Pursuit is initialized with 'guess ROIs', which it then refines. You need to pass these guess ROIs to spikePursuit.m as the third argument.
To see how it works I recommend loading the ROI file in the exampledata folder (using the loadROIs scripts) and calling
spikePursuit([],[], ROIs) then selecting the data block in the exampledata folder. 


