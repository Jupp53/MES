# MES -  METODA ELEMENTÓW SKOŃCZONYCH
# FEM -  FINITE ELEMENT METHOD
The task of the program was to conduct a simulation of the unsteady heat flow process with a boundary condition of convection. To achieve this, Fourier's equation for an unsteady process was used.

The entire project was written in C++. The majority of the program is written in the main function, step by step leading to the solution of the studied problem. The program contains structures, namely:

* struct GlobalData - Contains data about the conditions for conducting the simulation obtained from a file.
* struct Node - A structure containing the coordinates of a node and its BC flag, indicating whether a boundary condition exists on its surface.
* struct Element - Has an array with the IDs of its nodes.
* struct Side - The structure of an element's side. It contains an array of integration points, ksi and eta coordinates, and an array of shape function values. In its constructor, the integration point coordinates on the element surface and the shape function values for each integration point of that side are calculated.
* struct Grid - The structure of the mesh. It contains an array of nodes and elements that make up the mesh and their number.

During data reading from the file, individual elements and their constituent nodes are created.
Then, for each element of the mesh, the following are calculated:
* tab_jacobians - Contains the already inverted Jacobian matrix. The Jacobian matrix is responsible for carrying out the geometric transformation, which allows the element to be transformed from the local system to the global system.
* det - The determinant of detJ is the ratio of the length of the global side to the length of the local side of the element. Used in the calculations of the H and C matrices.
* matrix_H - The matrix of finite element characteristics. It shows how heat transfer occurs in the element.
* matrix_C - Indicates how much heat one finite element can store. It depends on the heat capacity and density.
In the C matrix, we use shape functions, and in the H matrix, we use their derivatives. Integration over dV is carried out by multiplying the result by the Jacobian of the transformation of that integration point.

* det_BC - The determinant of the boundary condition Jacobian equals L/2, where L is the length of the wall segment (calculated using the Pythagorean theorem), and 2 is the length of the segment in the local system.
* matrixHBC - Takes into account the boundary condition, i.e., it describes the behavior of the material at the boundary. In this case, the integral is calculated over the surface, just like in the P vector. In the case of the created program, the condition is calculated for each side, and then, if there is no boundary condition on it, the matrix is set to zero.
* vector_P - Describes the intensity of external forces acting on the system.
Later, the HBC matrix is added to the H matrix of individual elements, and then the aggregation of the H, C matrices, and P vector is started.

After this stage, the system of equations is solved using the Gaussian elimination method.
* bool gauss(int n, double** AB, double* X) - The Gaussian elimination method
* MinTemp and MaxTemp - Arrays containing the minimum and maximum values of the sought temperatures at the nodes after each iteration.
