#include "pch.h"
#include "Stitch.h"
#define DEBUG_NICP 1

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::SparseMatrix;

/* IO part */

void printSparse(SparseMatrix<double> M,std::fstream &stream)
{
	int p = M.nonZeros();
	stream << "Printing Matrix" << std::endl;
	stream << M.nonZeros() << std::endl;
	stream << M.rows() << "  " << M.cols() << std::endl;
	for (Eigen::Index k = 0; k < M.outerSize(); k++)
	{
		for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
		{
			
			stream << it.row() << "  " << it.col() << "  "<< it.value()<< std::endl;
		}
	}
}

void printSparse(SparseMatrix<double, Eigen::RowMajor> M, std::fstream& stream)
{
	int p = M.nonZeros();
	stream << "Printing Matrix" << std::endl;
	stream << M.nonZeros() << std::endl;
	stream << M.rows() << "  " << M.cols() << std::endl;
	for (Eigen::Index k = 0; k < M.outerSize(); k++)
	{
		for (SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(M, k); it; ++it)
		{
			stream << it.row() << "  " << it.col() << "  " << it.value() << std::endl;
		}
	}
}

Stitch::FileReader::FileReader(std::string filename):data()
{
	std::fstream fS(filename);
	if (fS.fail())
	{
		std::cout << filename << std::endl;
		throw;
	}
	int a, b;
	fS >> a >> b;
	MatrixXd m(a, b);
	for (int i = 0; i < a; i++)
	{
		for (int j = 0; j < b; j++)
		{
			double x;
			fS >> x;
			m(i, j) = x;
		}
	}
	data = m;
	fS.close();
	

}

MatrixXd Stitch::FileReader::fread(std::string filename)
{
	std::fstream fS(filename);
	int a, b;
	if (fS.fail())
	{
		std::cout << filename << std::endl;
		throw;
	}
	fS >> a >> b;
	MatrixXd m(a, b);
	for (int i = 0; i < a; i++)
	{
		for (int j = 0; j < b; j++)
		{
			double x;
			fS >> x;
			m(i, j) = x;
		}
	}
	data = m;
	fS.close();
	return data;
}

Stitch::ObjReader::ObjReader(std::string filename)
{
	using std::vector;
	std::fstream fS(filename);
	char c;
	// Vectors for saving postion
	vector<double> vX;
	vector<double> vY;
	vector<double> vZ;

	struct face {
		face(int a, int b, int c) :a(a - 1), b(b - 1), c(c - 1) {};
		int a, b, c;
	};
	vector<face> vf;

	while (fS >> c)
	{
		double x, y, z;
		if (c == 'v')
		{
			fS >> x >> y >> z;
			vX.push_back(x);
			vY.push_back(y);
			vZ.push_back(z);
		}
		if (c == 'f')
		{
			int x, y, z;
			fS >> x >> y >> z;
			vf.push_back(face(x, y, z));

		}
	}
	data.resize(vX.size(), 3);
	for (int i = 0; i < vX.size(); i++)
	{
		data(i, 0) = vX[i];
		data(i, 1) = vY[i];
		data(i, 2) = vZ[i];
	}

	MatrixXi facedata;
	facedata.resize(vf.size(), 3);

	for (int i = 0; i < vf.size(); i++)
	{
		facedata(i, 0) = vf[i].a;
		facedata(i, 1) = vf[i].b;
		facedata(i, 2) = vf[i].c;
	}

	objmesh.vertices = data;
	objmesh.edges = facedata;
}


Stitch::RegressionHead::RegressionHead(MatrixXd head_mean, MatrixXd head_U, MatrixXd Whf, MatrixXd input_face, MatrixXd face_U, MatrixXd face_mean)
{
	std::cout << face_mean.rows() << " " << face_mean.cols() << std::endl;
	input_face.resize(input_face.rows() * input_face.cols(), 1);
	auto step1 = input_face - face_mean;
	auto step2 = face_U.transpose() * step1;
	auto step3 = Whf * step2;
	auto step4 = head_U * step3;
	auto result = head_mean + step4;
	data = result;
	
	return ;
}

/*  Mesh  */

Stitch::Mesh::Mesh(double v[], int vertices_len, int t[], int triangles_len)
{
	edges = MatrixXi(triangles_len / 3, 3);
	vertices = MatrixXd(vertices_len / 3, 3);
	for (int i = 0; i < vertices_len/3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			vertices(i, j) = v[i * 3 + j];
		}
	}
	for (int i = 0; i < triangles_len / 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			edges(i, j) = t[i * 3 + j];
		}
	}
}

Stitch::Mesh::Mesh(::std::vector<Eigen::Vector3f> V, ::std::vector<::std::array<int, 3>> F)
{
	vertices = MatrixXd(V.size(), 3);
	edges = MatrixXi(F.size(), 3);
	for (int i = 0; i < V.size(); i++)
	{
		vertices(i, 0) = V[i](0);
		vertices(i, 1) = V[i](1);
		vertices(i, 2) = V[i](2);
	}
	for (int i = 0; i < F.size(); i++)
	{
		edges(i, 0) = F[i][0];
		edges(i, 1) = F[i][1];
		edges(i, 2) = F[i][2];
	}
}

Stitch::mesh_v_and_f Stitch::Mesh::ToEos()
{
	mesh_v_and_f result;
	for (int i = 0; i < vertices.rows(); ++i) {
		result.vertices.push_back(Eigen::Vector3f(vertices(i, 0), vertices(i, 1), vertices(i, 2)));
	}
	for (int i = 0; i < edges.rows(); ++i)
	{
		result.triangles.push_back(std::array<int, 3>({ edges(i,0),edges(i,1),edges(i,2) }));
	}
	return result;
}

void Stitch::Mesh::export_as_obj_file(std::string filename)
{
	std::fstream filestream(filename, ::std::ios::out);
	if (!filestream.is_open())
	{
		std::cout << "file not open!" << std::endl;
	}
	for (int i = 0; i < vertices.rows(); i++)
	{
		filestream << "v " << vertices(i, 0) << ' ' << vertices(i, 1) << ' ' << vertices(i, 2) << std::endl;
	}

	for (int j = 0; j < edges.rows(); j++)
	{
		filestream << "f " << edges(j, 0) + 1<< ' ' << edges(j, 1) + 1 << ' ' << edges(j, 2) + 1 << std::endl;
	}
	filestream.close();
	return;
}

void Stitch::Mesh::export_boudary_out(std::string filename, MatrixXi Ref)
{
	std::fstream filestream(filename, 'w');
	if (filestream.is_open())
	{
		for (int i = 0; i < Ref.rows(); i++)
		{
			vertices(Ref(i, 0), 2) += 100;
		}
		for (int i = 0; i < vertices.rows(); i++)
		{
			filestream << "v " << vertices(i, 0) << ' ' << vertices(i, 1) << ' ' << vertices(i, 2) << std::endl;
		}
		for (int j = 0; j < edges.rows(); j++)
		{
			filestream << "f " << edges(j, 0) + 1 << ' ' << edges(j, 1) + 1 << ' ' << edges(j, 2) + 1 << std::endl;
		}

	}
	else {
		std::cout << "File not open!!" << std::endl;
	}
	filestream.close();
	return;
}

void Stitch::Mesh::export_boudary_as_obj_file(std::string filename, MatrixXi Ref)
{
	std::fstream filestream(filename, 'w');
	if (filestream.is_open())
	{
		for (int k = 0; k < Ref.rows(); k++)
		{
			int i = Ref(k, 0);
			if (i == -1) continue;
			filestream << "v " << vertices(i, 0) << ' ' << vertices(i, 1) << ' ' << vertices(i, 2) << std::endl;
		}

	}
	else {
		std::cout << "File not open!!" << std::endl;
	}
	filestream.close();
	return;

}


void Stitch::Mesh::export_I(std::string filename, MatrixXi Ref)
{
	// Find a vertice on head
	int start_index = -1;
	for (int i = 0; i < vertices.rows(); i++)
	{
		if (vertices(i, 2) < -100)
		{
			start_index = i;
			break;
		}
	}
	if (start_index == -1)
	{
		std::cout << vertices.rows() << std::endl;
		std::cout << "Can not find a head vertices" << std::endl;
	}
	
	std::vector<std::set<int>> Adj(vertices.rows(), std::set<int>());
	// Make adjacent Matrix
	for (int i = 0; i < edges.rows(); i++)
	{
		int a = edges(i, 0);
		int b = edges(i, 1);
		int c = edges(i, 2);
		Adj[a].insert(b);
		Adj[a].insert(c);
		Adj[b].insert(a);
		Adj[b].insert(c);
		Adj[c].insert(a);
		Adj[c].insert(b);
	}
	// Add
	std::vector<int> back;
	back.push_back(start_index);
	std::vector<bool> searched(vertices.rows(), false);
	// Block
	for (int i = 0; i < Ref.rows(); i++)
	{
		if (Ref(i, 0) == -1) continue;
		searched[Ref(i, 0)] = true;
	}
	// Block2
	std::set<int> block2;
	for (int i = 0; i < edges.rows(); i++)
	{
		int a = edges(i, 0);
		int b = edges(i, 1);
		int c = edges(i, 2);
		int count = 0;
		count += (searched[a]) ? 1 : 0 + (searched[b]) ? 1 : 0 + (searched[c]) ? 1 : 0;
		if (count >= 1)
		{
			block2.insert(a);
			block2.insert(b);
			block2.insert(c);
		}
	}
	for (int i : block2)
	{
		searched[i] = true;
	}
	int x = 0;
	while (x < back.size())
	{
		bool changed = false;
		int node = back[x];
		if (!searched[node]) {
			for (auto con : Adj[node])
			{
				back.push_back(con);
				//searched[con] = true;
			}
		}	
		searched[node] = true;
		x++;
	}
	std::set<int> exloop;
	for (auto i : back)
	{
		for (auto j : Adj[i]) 
		{
			exloop.insert(j);
		}
	}
	for (auto i : exloop)
	{
		back.push_back(i);
	}
	// Add the vertices use queue
	std::vector<bool> backed(vertices.rows(), false);

	for (auto i : back)
	{
		backed[i] = true;
	}
	std::vector<int> col(vertices.rows(), 0);
	{
		int u = 0;
		for (int i = 0; i < vertices.rows(); i++)
		{
			col[i] = u;
			if (backed[i]) u++;
		}
	}
	//
	std::fstream fs(filename,'w');
	for (int i = 0; i < vertices.rows(); i++)
	{	
		if (backed[i]) {
			fs << "v " << vertices(i, 0) << ' ' << vertices(i, 1) << ' ' << vertices(i, 2) << std::endl;
		}
	}
	for (int i = 0; i < edges.rows(); i++)
	{
		if (backed[edges(i, 0)] && backed[edges(i, 1)] && backed[edges(i, 2)]) {
			fs << "f " << col[edges(i, 0)] + 1 << ' ' << col[edges(i, 1)] + 1 << ' ' << col[edges(i, 2)] + 1 << std::endl;
		}
		
	}
	fs.close();
}

void Stitch::Mesh::export_out(double*& v, int& vl, int*& t, int& tl)
{
	vl = vertices.rows() * vertices.cols();
	tl = edges.rows() * edges.cols();
	for (int i = 0; i < vertices.rows(); i++)
	{
		for (int j = 0; j < vertices.cols(); j++)
		{
			v[i * 3 + j] = vertices(i, j);
		}
	}
	for (int i = 0; i < edges.rows(); i++)
	{
		for (int j = 0; j < edges.cols(); j++)
		{
			t[i * 3 + j] = edges(i, j);
		}
	}
	return;
}

void Stitch::Mesh::M(MatrixXd facevertices, MatrixXi faceedges, MatrixXi ref)
{
	for (int i = 0; i < contour.size(); i++)
	{
		int x = ref(contour[i] - 1, 0);
		vertices(x, 0) = facevertices(contour[i] - 1, 0);
		vertices(x, 1) = facevertices(contour[i] - 1, 1);
		vertices(x, 2) = facevertices(contour[i] - 1, 2);
	}
}

void Stitch::Mesh::leave_back(MatrixXi Ref)
{
	// Find a vertice on head
	int start_index = -1;
	for (int i = 0; i < vertices.rows(); i++)
	{
		if (vertices(i, 2) < -100)
		{
			start_index = i;
			break;
		}
	}
	if (start_index == -1)
	{
		std::cout << "Can not find a head vertices" << std::endl;
	}

	std::vector<std::set<int>> Adj(vertices.rows(), std::set<int>());
	for (int i = 0; i < edges.rows(); i++)
	{
		int a = edges(i, 0);
		int b = edges(i, 1);
		int c = edges(i, 2);
		Adj[a].insert(b);
		Adj[a].insert(c);
		Adj[b].insert(a);
		Adj[b].insert(c);
		Adj[c].insert(a);
		Adj[c].insert(b);
	}
	// Add
	std::vector<int> back;
	back.push_back(start_index);
	std::vector<bool> searched(vertices.rows(), false);
	// Block
	for (int i = 0; i < Ref.rows(); i++)
	{
		if (Ref(i, 0) == -1) continue;
		searched[Ref(i, 0)] = true;
	}
	// Block2
	std::set<int> block2;
	for (int i = 0; i < edges.rows(); i++)
	{
		int a = edges(i, 0);
		int b = edges(i, 1);
		int c = edges(i, 2);
		int count = 0;
		count += (searched[a]) ? 1 : 0 + (searched[b]) ? 1 : 0 + (searched[c]) ? 1 : 0;
		if (count >= 1)
		{
			block2.insert(a);
			block2.insert(b);
			block2.insert(c);
		}
	}
	for (int i : block2)
	{
		searched[i] = true;
	}
	int x = 0;
	while (x < back.size())
	{
		bool changed = false;
		int node = back[x];
		if (!searched[node]) {
			for (auto con : Adj[node])
			{
				back.push_back(con);
			}
		}
		searched[node] = true;
		x++;
	}
	std::set<int> exloop;
	for (auto i : back)
	{
		for (auto j : Adj[i])
		{
			exloop.insert(j);
		}
	}
	for (auto i : exloop)
	{
		back.push_back(i);
	}
	// Add the vertices use queue
	std::vector<bool> backed(vertices.rows(), false);

	for (auto i : back)
	{
		backed[i] = true;
	}
	std::vector<int> col(vertices.rows(), 0);
	{
		int u = 0;
		for (int i = 0; i < vertices.rows(); i++)
		{
			col[i] = u;
			if (backed[i]) u++;
		}
	}

	int rows = 0;
	for (int i = 0; i < vertices.rows(); i++)
	{
		if (backed[i]) {
			rows++;
		}
	}
	int edges_rows = 0;
	for (int i = 0; i < edges.rows(); i++)
	{
		if (backed[edges(i, 0)] && backed[edges(i, 1)] && backed[edges(i, 2)]) {
			//fs << "f " << col[edges(i, 0)] + 1 << ' ' << col[edges(i, 1)] + 1 << ' ' << col[edges(i, 2)] + 1 << std::endl;
			edges_rows++;
		}

	}
	MatrixXd newvertices(rows, 3);
	int vi = 0;
	std::fstream fs("HeadVerticesList.txt",::std::ios::out);
	for (int i = 0; i < vertices.rows(); i++)
	{
		if (backed[i]) {
			newvertices(vi, 0) = vertices(i, 0);
			newvertices(vi, 1) = vertices(i, 1);
			newvertices(vi, 2) = vertices(i, 2);
			vi++;
			fs << i << ' ';
		}
		
	}
	fs.close();
	int ei = 0;
	MatrixXi newedges(edges_rows, 3);
	for (int i = 0; i < edges.rows(); ++i)
	{
		if (backed[edges(i, 0)] && backed[edges(i, 1)] && backed[edges(i, 2)])
		{
			newedges(ei, 0) = col[edges(i, 0)];
			newedges(ei, 1) = col[edges(i, 1)];
			newedges(ei, 2) = col[edges(i, 2)];
			ei++;
		}
	}
	vertices = newvertices;
	edges = newedges;
	//fs.close();
}

Stitch::NICP::NICP(Mesh source, Mesh target)
{
	KnnSearch knn(target.vertices); // KNN
	// Building the Upper part of Matrix A

	auto kron_M_G = building_A_Upper(source.edges,source.vertices);

	SparseMatrix<double,Eigen::RowMajor> D(source.vertices.rows(),source.vertices.rows() * 4);
	D = building_A_Lower_init(source.vertices);
	SparseMatrix<double,Eigen::RowMajor> X = building_X_init(source.vertices.rows());
	MatrixXi matchIndices;
	// Iteration match starts
	for (int i = 0; i < 20; i++)
	{
		double alphas = 200 - 199 * (i / 20.0);

		SparseMatrix<double,Eigen::RowMajor> vertsTransformed = D * X;
		matchIndices = knn.Kneighbors(vertsTransformed,1,225.0);
		SparseMatrix<double, Eigen::RowMajor> B_lower
			= braket_row(target.vertices.sparseView(),matchIndices);//?

		SparseMatrix<double,Eigen::RowMajor> U0(kron_M_G.rows(), 3);
		SparseMatrix<double,Eigen::RowMajor> DwVec = braket_remove(D, matchIndices);
		
		auto ArowBase = StackRows(alphas * kron_M_G, DwVec);
		auto BrowBase = StackRows(U0, B_lower);
		
		auto A = RowBaseToColBase(ArowBase);
		auto B = RowBaseToColBase(BrowBase);
		auto Xrow  =  CholeskySolve(A, B);
		
		auto X2 = CholeskySolve(RowBaseToColBase(DwVec), RowBaseToColBase(B_lower));
		
		X = ColBaseToRowBase(Xrow);
	}
	SparseMatrix<double> results = D * X;// Result of NICP
	deformed_mesh.vertices = results;
	deformed_mesh.edges = source.edges;
	refindices = matchIndices;

	return;
}

void Stitch::NICP::get_relation_vertex_only(Mesh source, Mesh target)
{
	KnnSearch knn(target.vertices);
	refindices = knn.Kneighbors(source.vertices, 1, 1600.0);
}

void Stitch::NICP::solve(Mesh source, Mesh target)
{
	auto kron_M_G = building_A_Upper(source.edges, source.vertices);

	SparseMatrix<double, Eigen::RowMajor> D(source.vertices.rows(), source.vertices.rows() * 4);
	D = building_A_Lower_init(source.vertices);
	SparseMatrix<double, Eigen::RowMajor> X = building_X_init(source.vertices.rows());
	MatrixXi matchIndices = ReverseRef(refindices,source.vertices.rows());
	// Iteration match starts
	std::vector<double> alpha_list = {150,100,50,1,0.001};
	std::vector<std::string> filenameList = { "5","10","15","20","25" };
	int counti = 0;
	for (double i: alpha_list)
	{
		double alphas = i;

		SparseMatrix<double, Eigen::RowMajor> vertsTransformed = D * X;
		SparseMatrix<double, Eigen::RowMajor> B_lower
			= braket_row(target.vertices.sparseView(), matchIndices);//?

		SparseMatrix<double, Eigen::RowMajor> U0(kron_M_G.rows(), 3);
		SparseMatrix<double, Eigen::RowMajor> DwVec = braket_remove(D, matchIndices);

		auto ArowBase = StackRows(alphas * kron_M_G, DwVec);
		auto BrowBase = StackRows(U0, B_lower);
		auto A = RowBaseToColBase(ArowBase);
		auto B = RowBaseToColBase(BrowBase);
		auto Xrow = CholeskySolve(A, B);
		auto X2 = CholeskySolve(RowBaseToColBase(DwVec), RowBaseToColBase(B_lower));
		X = ColBaseToRowBase(Xrow);
			
		SparseMatrix<double> results = D * X;
		/*
		Mesh T;
		T.vertices = results;
		T.edges = source.edges;
		
		std::string testfilefolder = "C:/Users/Yueyuan/Documents/GitHub/nonrigid_icp/exp3/test1/";
		std::string testfilename = testfilefolder + filenameList[counti] + ".obj";
		std::string testfilename2 = testfilefolder + filenameList[counti] + "cc" + ".obj";
		counti++;
		T.M(target.vertices, target.edges, refindices);
		T.export_I(testfilename, refindices);
		T.export_as_obj_file(testfilename2);
		*/
	}

	SparseMatrix<double> results = D * X;// Result of NICP
	
	deformed_mesh.vertices = results.toDense();
	deformed_mesh.edges = source.edges;
	return;
}

void Stitch::NICP::solve(Mesh source, Mesh target, int step)
{
	if (step < 10) return;
	auto kron_M_G = building_A_Upper(source.edges, source.vertices);
	if (step < 11) return;
	SparseMatrix<double, Eigen::RowMajor> D(source.vertices.rows(), source.vertices.rows() * 4);
	if (step < 12) return;
	D = building_A_Lower_init(source.vertices);
	if (step < 13) return;
	SparseMatrix<double, Eigen::RowMajor> X = building_X_init(source.vertices.rows());
	if (step < 14) return;
	MatrixXi matchIndices = ReverseRef(refindices, source.vertices.rows());
	// Iteration match starts
	if (step < 15) return;
	std::vector<double> alpha_list = { 150,100,50,1,0.001 };
	if (step < 16) return;
	std::vector<std::string> filenameList = { "5","10","15","20","25" };
	if (step < 17) return;
	int counti = 0;
	//for (double i : alpha_list)
	{
		if (step < 18) return;
		double alphas = alpha_list[4];
		if (step < 19) return;
		SparseMatrix<double, Eigen::RowMajor> vertsTransformed = D * X;
		if (step < 20) return;
		SparseMatrix<double, Eigen::RowMajor> B_lower
			= braket_row(target.vertices.sparseView(), matchIndices);//?
		if (step < 21) return;
		SparseMatrix<double, Eigen::RowMajor> U0(kron_M_G.rows(), 3);
		if (step < 22) return;
		SparseMatrix<double, Eigen::RowMajor> DwVec = braket_remove(D, matchIndices);
		if (step < 23) return;
		auto ArowBase = StackRows(alphas * kron_M_G, DwVec);
		if (step < 24) return;
		auto BrowBase = StackRows(U0, B_lower);
		if (step < 25) return;
		auto A = RowBaseToColBase(ArowBase);
		if (step < 26) return;
		auto B = RowBaseToColBase(BrowBase);
		if (step < 27) return;
		auto Xrow = CholeskySolve(A, B);
		if (step < 28) return;
		auto X2 = CholeskySolve(RowBaseToColBase(DwVec), RowBaseToColBase(B_lower));
		if (step < 29) return;
		X = ColBaseToRowBase(Xrow);
		if (step < 30) return;
		SparseMatrix<double> results = D * X;
		/*
		Mesh T;
		T.vertices = results;
		T.edges = source.edges;

		std::string testfilefolder = "C:/Users/Yueyuan/Documents/GitHub/nonrigid_icp/exp3/test1/";
		std::string testfilename = testfilefolder + filenameList[counti] + ".obj";
		std::string testfilename2 = testfilefolder + filenameList[counti] + "cc" + ".obj";
		counti++;
		T.M(target.vertices, target.edges, refindices);
		T.export_I(testfilename, refindices);
		T.export_as_obj_file(testfilename2);
		*/
	}
	if (step < 31) return;
	SparseMatrix<double> results = D * X;// Result of NICP
	if (step < 32) return;
	deformed_mesh.vertices = results.toDense();
	deformed_mesh.edges = source.edges;
	if (step < 33) return;
	return;
}

MatrixXi Stitch::NICP::ReverseRef(MatrixXi mat,int len)
{

	MatrixXi result = -MatrixXi::Ones(len,1);
	for (int i = 0; i < mat.rows(); i++)
	{
		int k = mat(i, 0);
		if (contour_set.find(i + 1) == contour_set.end())
		{
			//continue;
		}
		if (k != -1) 
		{
			result(k, 0) = i;
		}
	}
	return result;
}

SparseMatrix<double,Eigen::RowMajor> Stitch::NICP::building_A_Upper(const MatrixXi &Edge, const MatrixXd &vertices)
{
	auto edgeSet = Stitch::GetEdgeSet(Edge);
	//std::cout << "Print edgeSet" << std::endl;
	for (auto U : edgeSet)
	{
		int a, b;
		std::tie(a, b) = U;
	}
	//std::cout << "Size of Edge Set " << edgeSet.size() << std::endl;
	auto M = EdgeMatrix(vertices, edgeSet);
	MatrixXd G = MatrixXd::Zero(4, 4);
	G(0, 0) = 1;
	G(1, 1) = 1;
	G(2, 2) = 1;
	G(3, 3) = 1;
	auto kron_M_G = Stitch::KroneckerProductDiagonal(M, G);

	return kron_M_G;
}

SparseMatrix<double,Eigen::RowMajor> Stitch::NICP::building_A_Lower_init(const MatrixXd& vertices)
{
	SparseMatrix<double, Eigen::RowMajor> D(vertices.rows(),vertices.rows() * 4);
	std::vector<Eigen::Triplet<double>> tripleList;
	for (int i = 0; i < vertices.rows(); i++)
	{
		tripleList.push_back(Eigen::Triplet<double>(i, i * 4,vertices(i, 0)));
		tripleList.push_back(Eigen::Triplet<double>(i, i * 4 + 1, vertices(i, 1)));
		tripleList.push_back(Eigen::Triplet<double>(i, i * 4 + 2, vertices(i, 2)));
		tripleList.push_back(Eigen::Triplet<double>(i, i * 4 + 3, 1)); // Correspende To Gamma		
	}
	D.setFromTriplets(tripleList.begin(), tripleList.end());
	return D;
}

SparseMatrix<double,Eigen::RowMajor> Stitch::NICP::building_X_init(int vertices_row)
{
	MatrixXd X_(4, 3);
	X_ << 1, 0, 0,
		0, 1, 0,
		0, 0, 1,
		0, 0, 0;
	SparseMatrix<double> X(4 * vertices_row, 3);
	std::vector<Eigen::Triplet<double>> tripleListX;

	for (int i = 0; i < vertices_row; i++)
	{
		tripleListX.push_back(Eigen::Triplet<double>(i * 4, 0, 1));
		tripleListX.push_back(Eigen::Triplet<double>(i * 4 + 1, 1, 1));
		tripleListX.push_back(Eigen::Triplet<double>(i * 4 + 2, 2, 1));
	}
	X.setFromTriplets(tripleListX.begin(), tripleListX.end());
	return X;
}

void Stitch::NICP::ExportRelation(std::string filename)
{
	std::fstream fs(filename,::std::ios::out);
	for (int i = 0; i < refindices.rows(); ++i)
	{
		fs << refindices(i, 0) << " ";
	}
	fs.close();
}

void Stitch::NICP::ImportRelation(::std::string filename)
{
	std::fstream fs(filename, ::std::ios::in);
	int x;
	std::vector<int> vs;
	while (fs >> x)
	{
		vs.push_back(x);
	}
	refindices = MatrixXi(vs.size(), 1);
	for (int i = 0; i < vs.size(); ++i)
	{
		refindices(i, 0) = vs[i];
	}
	
	fs.close();
}

SparseMatrix<double> Stitch::KroneckerProduct(SparseMatrix<double,Eigen::RowMajor> A, Eigen::MatrixXi B)
{	
	return SparseMatrix<double>();
}

SparseMatrix<double> Stitch::KroneckerProductDiagonal(SparseMatrix<double,Eigen::RowMajor> A, Eigen::MatrixXd B)
{
	SparseMatrix<double> result(A.rows() * B.rows(), A.cols() * B.cols());
	//std::cout << "Non Zeros of A " << A.nonZeros() << std::endl;
	std::vector<Eigen::Triplet<double>> T;
	
	for (int k = 0; k < A.outerSize(); ++k)
	{
		for (SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(A, k); it; ++it)
		{		
			for (int i = 0; i < B.rows(); i++)
			{
				if (B(i, i) < 0.001)
				{
					continue;
				}
				T.push_back(Eigen::Triplet<double>(it.row() * B.rows() + i, it.col() * B.cols() + i, it.value()));
			}
		}
	}
	result.setFromTriplets(T.begin(), T.end());
	//std::cout << "Non Zeros of result" << result.nonZeros() << std::endl;
	return result;
}

SparseMatrix<double, Eigen::RowMajor> Stitch::StackRows(SparseMatrix<double, Eigen::RowMajor> A, SparseMatrix<double, Eigen::RowMajor> B)
{
	if (A.cols() != B.cols())
	{
		std::cerr << "test" << std::endl;
		std::cerr << "Matrices with different cols can not stack together" << std::endl;
		std::cerr << A.cols() << std::endl;
		std::cerr << B.cols() << std::endl;
		throw;
		return SparseMatrix<double>();
	}
	SparseMatrix<double, Eigen::RowMajor> result(A.rows() + B.rows(), A.cols());
	result.reserve(A.nonZeros() + B.nonZeros());
	if (A.IsRowMajor && B.IsRowMajor)
	{
		for (Eigen::Index i = 0; i < A.outerSize(); i++)
		{
			result.startVec(i);
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator itA(A, i); itA; ++itA) {
				result.insertBack(itA.row(), itA.col()) = itA.value();
			}
		}
		for (Eigen::Index i = 0; i < B.outerSize(); i++)
		{
			result.startVec(i + A.outerSize());
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator itB(B, i); itB; ++itB)
			{
				result.insertBack(itB.row() + A.rows(), itB.col()) = itB.value();
			}
		}
	}
	//std::cout << "Stack rows successfully!! " << std::endl;
	if ((!A.IsRowMajor) && (!B.IsRowMajor))
	{
		for (Eigen::Index i = 0; i < A.outerSize(); i++)
		{
			result.startVec(i);
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator itA(A, i); itA; ++itA)
				result.insertBack(itA.row(), itA.col()) = itA.value();
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator itB(B, i); itB; ++itB)
				result.insertBack(itB.row() + A.rows(), itB.col()) = itB.value();
		}
	}
	result.finalize();
	//std::cout << "Stack Finalize Sucessfully!! " << std::endl;
	return result;
	
}

::std::set<::std::tuple<int, int>> Stitch::GetEdgeSet(MatrixXi faces)
{
	::std::set<::std::tuple<int, int>> result;
	for (int i = 0; i < faces.rows(); i++)
	{
		int a, b, c;
		a = faces(i, 0);
		b = faces(i, 1);
		c = faces(i, 2);
		if (b > a) {
			int t = a;
			a = b;
			b = t;
		}
		if (c > a)
		{
			int t = a;
			a = c;
			c = t;
		}
		if (c > b)
		{
			int t = b;
			b = c;
			c = t;
		}
		result.insert(::std::tuple<int, int>(a,b));
		result.insert(::std::tuple<int, int>(b,c));
		result.insert(::std::tuple<int, int>(a,c));
	}
	return result;
}

SparseMatrix<double,Eigen::RowMajor> Stitch::EdgeMatrix(MatrixXd vertices, EdgeSet edges)
{
	SparseMatrix<double> result(edges.size(),vertices.rows());
	std::vector<Eigen::Triplet<double>> T;
	int i = 0;
	for (auto edgeI : edges)
	{
		int a, b;
		std::tie(a, b) = edgeI;
		T.push_back(Eigen::Triplet<double>(i, a, -1));
		T.push_back(Eigen::Triplet<double>(i, b, 1));
		//result.insert(i, a) = -1;
		//result.insert(i, b) = 1;
		i++;
	}
	result.setFromTriplets(T.begin(), T.end());
	return result;
}

SparseMatrix<double,Eigen::RowMajor> Stitch::MulCofWise(SparseMatrix<double,Eigen::RowMajor> D, MatrixXd wVec)
{
	std::vector<Eigen::Triplet<double>> T;
	if (D.rows() == wVec.rows())
	{
		SparseMatrix<double> result(D.rows(), D.cols());
		for (int k = 0; k < D.outerSize(); ++k)
		{
			for (SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(D, k); it; ++it)
			{
				T.push_back(Eigen::Triplet<double>(it.row(),it.col(), it.value() * wVec(it.row(), 0)));
				//result.insert(it.row(), it.col()) =  it.value() * wVec(it.row(),0);
			}
		}
		result.setFromTriplets(T.begin(), T.end());
		return result;
	}
	else
	{
		std::cerr << "invalid rows number" << std::endl;
		return SparseMatrix<double>();
	}
}

SparseMatrix<double> Stitch::CholeskySolve(SparseMatrix<double> A, SparseMatrix<double> B)
{
	auto A_t = A.transpose();
	auto AtA = A_t * A;
	Eigen::SimplicialCholesky<SparseMatrix<double>> MatricesCholesky(AtA);
	auto XRHS = A_t * B;
	auto X = MatricesCholesky.solve(XRHS);
	return X;
}

SparseMatrix<double> Stitch::QRSolve(SparseMatrix<double> A, SparseMatrix<double>B) {
	Eigen::SparseQR<SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	solver.compute(A);
	SparseMatrix<double> ret = solver.solve(B);
	return ret;
}

SparseMatrix<double, Eigen::RowMajor> Stitch::braket_row(const SparseMatrix<double, Eigen::RowMajor> &X, const MatrixXi &ref)
{
	SparseMatrix<double,Eigen::RowMajor> result(ref.rows(), X.cols());
	std::vector<Eigen::Triplet<double>> tripleList;
	for (int i = 0; i < ref.rows(); i++)
	{
		int k = ref(i, 0);
		if (k == -1) continue;
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(X, k); it; ++it)
		{
			tripleList.push_back(Eigen::Triplet<double>(i, it.col(), it.value()));
		}
	}		
	result.setFromTriplets(tripleList.begin(), tripleList.end());
	return result;
}

SparseMatrix<double, Eigen::RowMajor> Stitch::braket_remove(const SparseMatrix<double, Eigen::RowMajor>& X, const MatrixXi& ref)
{
	SparseMatrix<double, Eigen::RowMajor> result(X.rows(), X.cols());
	std::vector<Eigen::Triplet<double>> tripleList;
	if (X.rows() == ref.rows())
	{
		for (int i = 0; i < ref.rows(); i++)
		{
			int k = ref(i, 0);
			if (k == -1) continue;
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(X, i); it; ++it)
			{
				tripleList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
			}
		}
	}
	else {
		std::cerr << "X rows() != ref.rows()" << std::endl;
		throw;
	}
	result.setFromTriplets(tripleList.begin(), tripleList.end());
	return result;
}

SparseMatrix<double> Stitch::RowBaseToColBase(SparseMatrix<double, Eigen::RowMajor> X)
{
	std::vector<Eigen::Triplet<double>> T;
	for (int i = 0; i < X.outerSize(); i++)
	{
		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(X, i); it; ++it)
		{
			T.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}
	SparseMatrix<double> result(X.rows(), X.cols());
	result.setFromTriplets(T.begin(), T.end());
	return result;
}

SparseMatrix<double, Eigen::RowMajor> Stitch::ColBaseToRowBase(SparseMatrix<double> X)
{
	std::vector<Eigen::Triplet <double>> T;
	for (int i = 0; i < X.outerSize(); i++)
	{
		for (SparseMatrix<double>::InnerIterator it(X, i); it; ++it)
		{
			T.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}
	SparseMatrix<double,Eigen::RowMajor> result(X.rows(), X.cols());
	result.setFromTriplets(T.begin(), T.end());
	return result;
}


Stitch::mesh_v_and_f Stitch::get__head(::std::vector<Eigen::Vector3f> vertice, ::std::vector<::std::array<int, 3>> triangles)
{
	//mesh_v_and_f result;
	Stitch::Mesh m(vertice,triangles);
	// Load model
	std::cout << m.vertices.rows() << " "  << m.vertices.cols();
	auto s = Fit(m);
	return s.ToEos();
}

Stitch::Mesh Stitch::Fit(Mesh m)
{
	// Load model
	std::string model_path = "";
	std::string folder = model_path;
	Stitch::FileReader Freader1(folder + "head_U.data");
	auto head_U = Freader1.data;
	auto head_mean = Freader1.fread(folder + "head_mean.data");
	auto Whf = Freader1.fread(folder + "Whf.data");
	auto face_U = Freader1.fread(folder + "face_U.data");
	auto face_mean = Freader1.fread(folder + "face_mean.data");

	// Compute the head
	Stitch::RegressionHead combiner(head_mean, head_U, Whf, m.vertices, face_U, face_mean);
	auto headdata = combiner.data;
	headdata.resize(3, headdata.size() / 3);
	headdata.transposeInPlace();
	
	
	std::string targetfile = folder + "Headv1.obj";
	Stitch::ObjReader targetreader(targetfile);
	auto x = targetreader.objmesh.edges;

	Stitch::Mesh m2;
	m2.vertices = headdata;
	m2.edges = x;

	
	//Stitch::NICP solver;
	//solver.ImportRelation();
	//solver.solve(m2, m);// fit m2(head) to m(face)
	//auto s = solver.deformed_mesh;
	//s.leave_back(solver.refindices);
	//auto re = merge(m, s);
	auto re = Stitch::Fit(m, m2);

	return re;
}

Stitch::Mesh Stitch::Fit(Mesh face, Mesh head)
{
	Stitch::NICP solver;
	solver.ImportRelation();
	solver.solve(head, face);
	auto s = solver.deformed_mesh;
	s.leave_back(solver.refindices);
	auto re = merge(face, s);
	return re;
}

Stitch::KnnSearch::KnnSearch(MatrixXd vertices)
{
	pointVec plist;
	for (int i = 0; i < vertices.rows(); ++i)
	{
		point_t point;
		for (int j = 0; j < vertices.cols(); ++j)
		{
			point.push_back(vertices(i, j));
		}
		plist.push_back(point);
	}
	kdtree = KDTree(plist);
}

MatrixXi Stitch::KnnSearch::Kneighbors(MatrixXd vertices,int K,double threhold)
{
	MatrixXi result(vertices.rows(), 1);
	for (int i = 0; i < vertices.rows(); i++)
	{
		point_t point;
		for (int j = 0; j < vertices.cols(); j++)
		{
			point.push_back(vertices(i, j));
		}
		auto nearest_index = kdtree.nearest_pointIndex(point);
		double dis = 0;
		auto res = nearest_index.first;
		for (int j = 0; j < vertices.cols(); j++)
		{
			dis += (res[j] - vertices(i, j))* (res[j] - vertices(i, j));
		}
		
		point_t pointP = nearest_index.first;
		 
		if (dis > threhold)
		{
			result(i, 0) = -1;
		}
		else {
			result(i, 0) = nearest_index.second;
		}
			
	}
	return result;
}

Stitch::Mesh Stitch::merge(Stitch::Mesh m1, Stitch::Mesh m2)
{
	Stitch::Mesh result;
	result.vertices = MatrixXd::Zero(m1.vertices.rows() + m2.vertices.rows(), m1.vertices.cols());
	result.edges = MatrixXi::Zero(m1.edges.rows() + m2.edges.rows(), m1.edges.cols());
	//result.vertices = Stitch::StackRows(m1.vertices, m2.vertices);
	for (int i = 0; i < m1.vertices.rows(); ++i)
	{
		for (int j = 0; j < m1.vertices.cols(); ++j)
		{
			result.vertices(i, j) = m1.vertices(i, j);
		}
	}
	for (int i = 0; i < m2.vertices.rows(); ++i) 
	{
		for (int j = 0; j < m2.vertices.cols(); ++j)
		{
			result.vertices(i + m1.vertices.rows(), j) = m2.vertices(i, j);
		}
	}

	for (int i = 0; i < m1.edges.rows(); ++i)
	{
		for (int j = 0; j < m1.edges.cols(); ++j)
		{
			result.edges(i, j) = m1.edges(i, j);
		}
	}

	for (int i = 0; i < m2.edges.rows(); ++i)
	{
		for (int j = 0; j < m2.edges.cols(); ++j)
		{
			result.edges(i + m1.edges.rows(), j) = m2.edges(i, j) + m1.vertices.rows();
		}
	}
	return result;
}

void __C__mesh__test()
{
	Stitch::ObjReader Oreader("C:/Users/Yueyuan/Documents/Exp/input_face.obj");
	auto facedata = Oreader.data;
	// Load model
	std::string folder = "C:/Users/Yueyuan/Documents/combining3Dmorphablemodels-master/combining3Dmorphablemodels-master/Prediction/";
	Stitch::FileReader Freader1(folder + "head_U.data");
	auto head_U = Freader1.data;
	auto head_mean = Freader1.fread(folder + "head_mean.data");
	auto Whf = Freader1.fread(folder + "Whf.data");
	auto face_U = Freader1.fread(folder + "face_U.data");
	auto face_mean = Freader1.fread(folder + "face_mean.data");

	Stitch::RegressionHead combiner(head_mean, head_U, Whf, facedata, face_U, face_mean);
	auto headdata = combiner.data;
	headdata.resize(3, headdata.size() / 3);
	headdata.transposeInPlace();

	std::string targetfile = "C:/Users/Yueyuan/Documents/GitHub/nonrigid_icp/exp3/Headv1.obj";
	Stitch::ObjReader targetreader(targetfile);
	auto x = targetreader.objmesh.edges;

	Stitch::Mesh m2;
	m2.vertices = headdata;
	m2.edges = x;

	//Set the solver 
	std::string sourcefile = "C:/Users/Yueyuan/Documents/GitHub/nonrigid_icp/exp3/Keeporder_pure.obj";
	//std::string targetfile = "C:/Users/Yueyuan/Documents/GitHub/nonrigid_icp/exp3/Headv1.obj";
	Stitch::ObjReader sourcereader(sourcefile);
	//Stitch::ObjReader targetreader(targetfile);

	Stitch::NICP solver(sourcereader.objmesh, targetreader.objmesh);
	solver.ExportRelation();
	solver.solve(m2,Oreader.objmesh);
	solver.deformed_mesh.export_as_obj_file("C:/Users/Yueyuan/Documents/Exp/K.obj");
	solver.deformed_mesh.export_I("C:/Users/Yueyuan/Documents/Exp/K_2.obj",solver.refindices);



	Stitch::Mesh mface = Oreader.objmesh;
	Stitch::ObjReader Oreader2("C:/Users/Yueyuan/Documents/Exp/K_2.obj");
	auto re = merge(mface, Oreader2.objmesh);
	re.export_as_obj_file("C:/Users/Yueyuan/Documents/Exp/Output.obj");
	
	return;

}

void __C__mesh__test2()
{
	Stitch::ObjReader Oreader("C:/Users/Yueyuan/Documents/Exp/Input_face.obj");
	Stitch::Mesh m = Oreader.objmesh;
	std::cout << m.vertices.rows() << " " << m.vertices.cols();
	auto s = Fit(m);
	s.export_as_obj_file("C:/Users/Yueyuan/Documents/Exp/K_4.obj");
	
	return;
}

void neck_test()
{
	Stitch::ObjReader Oreader("C:/Users/Yueyuan/Documents/Exp/Input_face.obj");
	Stitch::Mesh m = Oreader.objmesh;
	std::cout << m.vertices.rows() << " " << m.vertices.cols();

	// Load model
	std::string model_path = "";
	std::string folder = model_path;
	Stitch::FileReader Freader1(folder + "head_U.data");
	auto head_U = Freader1.data;
	auto head_mean = Freader1.fread(folder + "head_mean.data");
	auto Whf = Freader1.fread(folder + "Whf.data");
	auto face_U = Freader1.fread(folder + "face_U.data");
	auto face_mean = Freader1.fread(folder + "face_mean.data");

	// Compute the head
	Stitch::RegressionHead combiner(head_mean, head_U, Whf, m.vertices, face_U, face_mean);
	auto headdata = combiner.data;
	headdata.resize(3, headdata.size() / 3);
	headdata.transposeInPlace();

	std::string targetfile = folder + "Headv1.obj";
	Stitch::ObjReader targetreader(targetfile);
	auto x = targetreader.objmesh.edges;

	Stitch::Mesh m2;
	m2.vertices = headdata;
	m2.edges = x;

	m2.export_as_obj_file("C:/Users/Yueyuan/Documents/Exp/HeadY.obj");
	//
	Stitch::ObjReader neckreader("C:/Users/Yueyuan/Documents/Exp/neckCirclePure.obj");
	

	Stitch::NICP solver;
	solver.get_relation_vertex_only(neckreader.objmesh, m2);
	solver.ExportRelation("NeckRelation.data");
	solver.solve(m2, neckreader.objmesh);// fit m2(head) to m(neck)
	auto s = solver.deformed_mesh;
	s.export_as_obj_file("FitNeck.obj");
	s.export_boudary_as_obj_file("NectBoundary.obj", solver.refindices);
	m2.export_boudary_as_obj_file("HeadNeck.obj", solver.refindices);

}