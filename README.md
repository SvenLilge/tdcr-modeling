## Modeling of Tendon Driven Continuum Robots

This code implements different approaches to model the kinematics/statics of a tendon driven continuum robot (TDCR) found in the current state of the art.
The implementation considers a two segment TDCR with three tendons per segment subject to an external load the robot's tip.
Both C++ and MATLAB code is provided.
The following modeling approaches are implemented:

- A geometry based modeling approach assuming constant curvature bending for each segment (Constant Curvature Model)
- A static modeling approach modeling assuming constant curvature bending for each subsegment, i.e. between neighbouring disks (Piecewise Constant Curvature Model)
- A static modeling approach modeling each subsegment as a pseudo rigid body with virtual discrete joints (Pseudo Rigid Body Model)
- A static modeling approach modeling each segment as a Cosserat rods subject to tendon loads (Cosserat Rod Model)
- A static modeling approach modeling each subsegment, i.e. between neighbouring disks, as a Cosserat rods subject to tendon loads (Subsegment Cosserat Rod Model)

Every model is based upon a state of the art modeling approach.
This repository is part of the following publication, which also provides theoretical derivations of each model as well as references to the original publications of each:

How to model tendon-driven continuum robots and benchmark modelling performance  
Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs  
frontiers in Robotics and AI 2021  
DOI: 10.3389/frobt.2020.630245  

Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

### Dependencies (C++)

The C++ implementation makes heavy use of both the Eigen library and the GNU Scientific Library, which need to be installed in order to compile the provided C++ code.

- [Eigen Library](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)
 
### Installation Instructions using CMake (C++)

In the root directory of the C++ code (where the main.cpp file is located) run the following commands:

	mkdir build
	cd build
	cmake ..
	make

Alternatively, the code can be compiled in Release mode for performance:
	
	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	make

Afterwards, you can execute the code by running the executable "tdcr-modeling".

### More Information

It might be worthwhile to dry different numerical solvers for the non-linear least squares fitting to check which offers the best convergence performance for a given set of TDCR parameters, tendon tensions and external loads. Different solvers can be commented and uncommented in the .cpp files of each implemented model. From experience, both the Levenberg-Marquadt solver (selected by default) and the Dogleg solver perform generally well.

If you found the provided implementation of the TDCR modeling approaches helpful or used parts of it yourself, please refer to it using the following BibTeX entry:

	@ARTICLE{RaoTDCRModelling,
	AUTHOR={Rao, Priyanka and Peyron, Quentin and Lilge, Sven and Burgner-Kahrs, Jessica},
	TITLE={How to Model Tendon-Driven Continuum Robots and Benchmark Modelling Performance},
	JOURNAL={Frontiers in Robotics and AI},
	VOLUME={7},
	YEAR={2021},
	DOI={10.3389/frobt.2020.630245},
	ISSN={2296-9144}
	}

