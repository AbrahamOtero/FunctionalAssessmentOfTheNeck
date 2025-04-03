First, the ExtractAngles.m script must be executed to extract the HMD positions for each rotaion axis from the raw login dat of the VR setup. This will generate new files under directories with the suffix _preprocessed.

Next, the SignalAnalysis.m script must be executed. This script calculates the spectrograms for the movement in each axis, and for the total movement. Then it applies statistical tests on the spectral power in each band, and on the duration of the records.

Optionally, the Spectrogram.m script can also be executed, which will generate .png files with the spectrograms of all subjects under the _preprocessed directory hierarchy.

The rest of the files in the repository are auxiliary functions for these three scripts.

Distributed under the MIT License (See accompanying file LICENSE or copy at http://opensource.org/licenses/MIT)
