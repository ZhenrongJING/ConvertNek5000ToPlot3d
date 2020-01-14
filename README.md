Nek5000 uses unstructured hexahedron mesh. But sometimes it is perferable to convert the mesh to structured one in post-processing. These codes are written for this purpose.

Attention: not all the unstructured mesh can be converted to structured one!

For now, the code only supports 3D elements. It also requires the user have a good knowlege on the distribution of elements. The usage of the code is demonstrated using an example flow around airfoil.

The above figure shows the simulation domain. The span of the airfoil is aligned with z direction. The span has 3 elements and periodic boundary condition is used at two span ends for the simulation. When is mesh is generated in ICEM, the domain is divided into 5 blocks, which are marked in the figure. 

The exported .msh file is converted to .rea file using mshconvert. After simulation finishes, we obtain the airfoil0.f000* files. 

1) compile and run read_mesh.f90 to read airfoil0.f0000*, this will output several .xyz file to mesh folder. Each .xyz file contains all the elements belongs to this blocks. (.xyz file is in plot3d format and can be opend with Tecplot.).

2) compile and run meshPressure.f90, which will output combined_pressure.xyz to combined_mesh folder. If we open combined_pressure.xyz file in Tecplot, it can be seen that now there is only one large structured mesh block. 

3) similar to reading mesh, we can use the same way to read velocity field and pressure.
$ gfortran read_uvw.f90 -o read_uvw
$ ./read_uvw 2 
