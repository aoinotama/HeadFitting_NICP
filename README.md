# HeadFitting_NICP
 
This the library code using NICP to Fit head to a certain face. 

## Configuration

Install C++ (https://eigen.tuxfamily.org/index.php?title=Main_Page) Eigen first.
Copy files in data/ folder to your working directory. And you can change the folder in the code too.
 
## API and function 
 
You can use the get__head function which input vertices and faces, it will return the vertices and faces of the Fitted Head.</br>

    mesh_v_and_f get__head(::std::vector<Eigen::Vector3f> vertice, ::std::vector<::std::array<int, 3>> triangles);
	Mesh Fit(Mesh x);

If you provide your own head, you can use 

	Mesh Fit(Mesh face, Mesh head);	

Use the constructor to construct the Mesh:

	Mesh(::std::vector<Eigen::Vector3f> vertices, ::std::vector<::std::array<int, 3>> triagles);

Or Use the Objreader to read the obj file(Not fully supported, please only keep the vertices lines and faces lines):

	class ObjReader
	{
	public:
		ObjReader(std::string filename);
		MatrixXd data;
		Mesh objmesh;
	};
	
The Mesh will be saved in .objmesh member.
