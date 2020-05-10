# SFND_Radar_Target_Generation_And_Detection
Udacity Sensor Fusion Radar Repository

#### Implementation of 2D CFAR
There's a total of 4 nested for-loops, the outer 2 for loops iterate through all the eligible rows and cols of RDM to select the CUT cell. While ignoring the outer perimeters to account for training and guard cells
The inner 2 loops will iterate through the box that includes the train, guard, and CUT cells to add the value of the RDM to the noise_level. An Ifelse statement is used to prevent guard and CUT cells to be included in the addition.
After iterating through the innter 2 loops, the summation of noise level is normalised by dividing it with number of training cells

#### Selection of Training, Guard cells and offset
The aim of the guard cells are to prevent the signal to leak to training cells, thus increasing the noise threshold excessively.
The raw RDM plot shows that the StDev of the velocity signal is a lot higher than range. Thus the training and guard cells for the velocity must be sufficiently high.
Since the range signal is narrow, a small training cells of 4 and guard cells of 2 is selected. While for velocity, a larger training cells of 8 and guard cells of 4 are chosen.
Although increasing the training and guard cells of the velocity will make the noise threshold better, it will reduce the range of velocity that we can detect due to the padding at the outer cols of the matrix

Offset is selected to clear the random noise generated, but not too high such that it is near the peaks of the real signal

#### Steps taken to suppress the non-thresholded cells at the edges
non-threshold cells at the edges are made to be zero by initialising the filteredRDM to be matrix of zeros with size(RDM)
