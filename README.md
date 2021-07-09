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
## Reference 
https://github.com/nabeel3133/combining3Dmorphablemodels
https://github.com/saikiran321/nonrigid_icp


## License
MIT License

Copyright (c) 2021 aoinotama

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
