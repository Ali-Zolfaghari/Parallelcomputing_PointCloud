# Parallelcomputing_PointCloud
Consider a point cloud that contain a number of spheres with definite coordinates. Naturally, each of these spheres will be surrounded by a number of points. First we want to identify the points that are in at least one of the spheres. On each point in space, a sphere of influence must be calculated as the value of a function (distance of two point in 3D coordinate). The value of this function depends on the coordinates of other points that are within its range. 
- A) First, an algorithm is written that obtains a function by dividing the work equal to or approximately equal to the number of processors. 
- B) The written parallel program reads the coordinate values of the spheres and points from two separate files, determines the points inside the sphere, calculates the radius of the impact range based on the fact that the minimum number of other points is in the impact limit of each point. Next, the function calculates for each point the average distance between the points in the range of impact for that point. Finally, the values of the function along with the coordinates of the points are stored in a file.
