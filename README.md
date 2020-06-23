# Maximum Inscribed Circle

## About
This is an OpenJUMP plugin designed to delimit Maximum Inscribed Circle (MIC) that can be inscribed in a polygon. It is based on the MICGIS algorithm developed by Beyhan et al. (2020). MICGIS algorithm performs an initial approximation of the medial axis (MA) of the polygon by using the Voronoi diagrams created for the points located along the polygon, and subsequently, makes use of the analytical/geometrical properties of the polygon, its edges and the vertices of approximate MA in conjunction with their general spatial configuration in order to delimit MIC of the polygon. In this process, MICGIS also benefits from the solutions proposed for the special cases of Apollonius’ Problem.

The first step in the algorithm is the removal of the points creating collinear segments along polygon edges. This is followed by the extraction of the points located along the polygon for the creation of input points for the construction of Voronoi diagram. During this process, extra points are created between the couple of points forming the edges of the polygon for a better approximation of the MA. Enrichment of the input points is followed by the creation of Voronoi diagrams for the input points. As the vertices of the Voronoi diagram inside the polygon provides us with an approximation of the MA, these vertices are considered as the candidates for the first approximation of MIC for the polygon in the algorithm.

In this respect, the first candidate for the center of MIC is revealed in the algorithm by maximizing the radius of circle tangent interior to one of the edges of polygon on the base of these vertices. The first candidate for the center of MIC is also used by the algorithm for the determination of the possible cores that may cover the center of MIC. For each core, MICGIS performs three approximation of MIC, and if necessary, the solutions proposed for the special cases of Apollonius’ Problem are also used for the approximation of the best MIC.

Beyhan, B., Güler, C. & Tağa, H. (2020) An algorithm for maximum inscribed circle based on Voronoi diagrams and geometrical properties. Journal of Geographical Systems, 22, 391–418. https://doi.org/10.1007/s10109-020-00325-3

## Installation and Use
In order to install “Maximum Inscribed Cirle” plugin, file named “MICGIS.jar” should be copied to the extension folder (“ext”) of OpenJUMP installation. In a typical Windows-based operating system, the default address for the respective folder is as following;

C:\Program Files\OpenJUMP-1.14.1-r6147-PLUS\lib\ext

Further instructions for the use and installation of the plugin can be found in "MICGIS User Guide.pdf" that can be downloaded from this website.

## Software Requirements
OpenJUMP - http://www.openjump.org/

## References

BEYHAN, B., GÜLER, C. & TAĞA, H. (2020) An algorithm for maximum inscribed circle based on Voronoi diagrams and geometrical properties. Journal of Geographical Systems, 22, 391–418. https://doi.org/10.1007/s10109-020-00325-3
