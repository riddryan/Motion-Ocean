# Motion-Ocean

Motion Ocean lets you create Inverse Dynamics models of humans from motion and force data in order to study the kinematics, dynamics, and energetics of human movements.  Unlike other Inverse Dynamics software, tracking of motion markers and the assignment of motion markers to specific body parts can be handled automatically.  Additionally, Motion Ocean can detect different types of locomotion and provide activity-specific analysis.  For example, if the data is of someone walking, it can detect when gait events such as heel-strikes and toe-offs, and segment the data according to those events in order to analyze the data on a per-step basis.

[TOC]

## Software ##

You need [Matlab 2014b](http://www.mathworks.com/products/matlab/) or later, with the global optimization toolbox.

## Getting Started ##
1. Clone the repository on your computer.
2. Download the contents of the folder **MotionOcean** from the HBCL server at [\\hbcl-server.engin.umich.edu\riddryan\MotionOcean](\\hbcl-server.engin.umich.edu\riddryan\MotionOcean).  Place these contents in the same folder that this (Running Perturbation) repository is in.
3. Clone the [HBCL toolbox repository](https://bitbucket.org/hbcl/hbcl-toolbox) on your computer.
4. Open Matlab and add the HBCL toolbox repository folder and it's subfolders to your matlab path.

## Classification of Human Models ##
Motion Ocean uses skeletal models from previous experiments to train an SVM (Support Vector Machine) classifier that maps markers to their location on the human body.  If you have new motion data, the  classifier will tell you if the marker is attached to the foot, the ankle, knee, pelvis, etc.  The classifier was trained on people assuming a standing pose with their arms outstreched.

### Feature Extraction ###
Motion capture data is loaded in a .c3d format and then processed before it is classified.  The pre-processing steps are:

1. Scale the position vectors of each marker such that their mean is equal to 1.
	* This should account for difference in height and width between subjects.  It may be better to scale the z-axis (vertical) separately.
2. Rotate the markers about the z-axis to minimize the covariance between the x &  y axes
	* This should take care of differences in the direction that subjects face during the standing trial

Each marker is a 6x1 feature vector where
1.  The first 3x1 features are equal to the 3D position of the marker relative to the center of position of all of the markers
2.  The last 3x1 features are equal to the ratio of the number of markers greater than the marker in each dimension compared to the number of total markers.


