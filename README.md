# Neuroscience Connectome Videos
Matlab Code to create Neuroscience Connectome videos.

The code is largely unorganised although all code needed to produce the videos are here.

There are 2 files here.

1) `sidevid.m`
2) `grabfiles.m`

The `sidevid.m` matlab file produces the video of the files side by side. It calls `grabfiles.m` and parses the relevant details for grabbing the information. Depending on the data being used this may be useless and it will be easier to read in the data directly into sidevid. 

The minimum data that sidevid requires is XYZ/MNI coordinates (labelled as `MNI` in `sidevid`), and time series data for each node with the corresponding voltage potentials (or whatever values you want to visualise).

The `e` values were used such that I had consistent structure across both images, for example if the left hand connectome was missing the 3rd node then it also makes sure to drop the 3rd node from the right hand connectome.

The information is also normalised as per the `if` statement on line 26. Set `normalise = False` to not normalise the data (i.e. mean = 0, sd = 1).

I warn that this code was written to be somewhat specific for the problem that was being considered at the time of writing, and a little bit of reworking would be required. I also have code purely for a single video but it requires a few touch ups.


