#pragma once

#ifdef STITCH_EXPORTS 
#define STITCH_API __declspec(dllexport)
#else
#define STITCH_API __declspec(dllimport)
#endif

#include<string>
#include<set>
#include<tuple>
#include<vector>
#include<array>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <deque>
#include <set>
#include <tuple>
#include "KDTree.hpp"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/SparseCore"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/SparseQR"
#include "KDTree.hpp"

#define __THREHOLD  225 


const std::vector<int> contour = { 12,22,23,24,26,39,41,357,359,360,361,364,365,366,370,371,372,374,381,382,383,385,386,387,389,392,393,394,442,455,465,466,467,476,478,774,776,777,778,781,782,783,787,788,789,791,796,797,798,800,801,802,803,806,807,808,1767,1772,1823,1827,1863,1864,1866,1867,1910,1913,1964,1967,2002,2005,2073,2076,2077,2081,2083,2086,2095,2096,2129,2131,2162,2164,2205,2208,2266,2269,2298,2301,2312,2313,2330,2333,2419,2421,2423,2424,2460,2463,2509,2510,2535,2536,2580,2581,2589,2590,2591,2592,2658,2659,2817,2818 };
const std::set<int> contour_set = { 12,22,23,24,26,39,41,357,359,360,361,364,365,366,370,371,372,374,381,382,383,385,386,387,389,392,393,394,442,455,465,466,467,476,478,774,776,777,778,781,782,783,787,788,789,791,796,797,798,800,801,802,803,806,807,808,1767,1772,1823,1827,1863,1864,1866,1867,1910,1913,1964,1967,2002,2005,2073,2076,2077,2081,2083,2086,2095,2096,2129,2131,2162,2164,2205,2208,2266,2269,2298,2301,2312,2313,2330,2333,2419,2421,2423,2424,2460,2463,2509,2510,2535,2536,2580,2581,2589,2590,2591,2592,2658,2659,2817,2818 };



namespace Stitch
{
	typedef ::std::set<::std::tuple<int, int>> EdgeSet;
	using Eigen::MatrixXd;
	using Eigen::MatrixXi;
	using Eigen::SparseMatrix;

	struct mesh_v_and_f {
		::std::vector<Eigen::Vector3f> vertices;
		::std::vector<::std::array<int, 3>> triangles;
	};

	/* Read matrix from files and save them into matrix
	*/
	class Mesh
	{
	public:
		MatrixXd vertices;
		MatrixXi edges;
		Mesh(double vertices[], int vertices_len, int triangles[], int triangles_len);
		Mesh() {};

		Mesh(::std::vector<Eigen::Vector3f> vertices, ::std::vector<::std::array<int, 3>> triagles);
		mesh_v_and_f ToEos();

		void export_as_obj_file(std::string filename);
		void export_boudary_out(std::string filename, MatrixXi Ref);
		void export_boudary_as_obj_file(std::string filename, MatrixXi Ref);
		void export_I(std::string filename, MatrixXi Ref);
		void export_out(double*& v, int& vl, int*& t, int& tl);
		void M(MatrixXd facevertices, MatrixXi faceedges, MatrixXi ref);
		void leave_back(MatrixXi Ref);
		
	};
	Mesh merge(Stitch::Mesh m1, Stitch::Mesh m2);
	class FileReader
	{
	public:
		FileReader(std::string filename);
		Eigen::MatrixXd data;
		Eigen::MatrixXd fread(std::string filename);
	};

	class ObjReader
	{
	public:
		ObjReader(std::string filename);
		MatrixXd data;
		Mesh objmesh;
	};

	class RegressionHead
	{
	public:
		RegressionHead(MatrixXd head_mean, MatrixXd head_U, MatrixXd Whf, MatrixXd input_face, MatrixXd faceU, MatrixXd face_mean);
		MatrixXd data;
	};

	class NICPresult {
	public:
		NICPresult(Mesh source1, Mesh target1);

		// 1. matches --> indices matches 
		// 2. deformed_mesh --> indices matches
	};

	class NICP {
	public:
		NICP(Mesh source, Mesh target);

		NICP() = default;
		void get_relation_vertex_only(Mesh source, Mesh target);
		
		void solve(Mesh source, Mesh target);
		Mesh resultMesh;
		Mesh deformed_mesh;
		MatrixXi refindices;
		::std::string model_folder;
		//Some const
		const int __nose__index = 114;
		const int __chin__index = 33;
		void solve(Mesh source, Mesh target, int step);
		MatrixXi ReverseRef(MatrixXi mat,int len);
		SparseMatrix<double,Eigen::RowMajor> building_A_Upper(const MatrixXi &Edge,const MatrixXd &vertices);
		SparseMatrix<double,Eigen::RowMajor> building_A_Lower_init(const MatrixXd& vertices);
		SparseMatrix<double, Eigen::RowMajor> building_X_init(int vertices_row);
		void ExportRelation(::std::string filename = "Relation.data");
		void ImportRelation(::std::string filename = "Relation.data");
	};

	SparseMatrix<double> KroneckerProduct(SparseMatrix<double,Eigen::RowMajor> A, Eigen::MatrixXi B);

	SparseMatrix<double> KroneckerProductDiagonal(SparseMatrix<double,Eigen::RowMajor> A, Eigen::MatrixXd B);

	SparseMatrix<double,Eigen::RowMajor> StackRows(SparseMatrix<double,Eigen::RowMajor> A, SparseMatrix<double,Eigen::RowMajor> B);

	::std::set<::std::tuple<int, int>> GetEdgeSet(MatrixXi faces);

	SparseMatrix<double, Eigen::RowMajor> EdgeMatrix(MatrixXd vertices, EdgeSet edges);

	SparseMatrix<double,Eigen::RowMajor> MulCofWise(SparseMatrix<double,Eigen::RowMajor> D, MatrixXd wVec);

	SparseMatrix<double> CholeskySolve(SparseMatrix<double> A, MatrixXd B);

	SparseMatrix<double> QRSolve(SparseMatrix<double> A, SparseMatrix<double> B);
	
	SparseMatrix<double> CholeskySolve(SparseMatrix<double> A, SparseMatrix<double> B);

	SparseMatrix<double, Eigen::RowMajor> braket_row(const SparseMatrix<double, Eigen::RowMajor> &X, const MatrixXi &ref);

	SparseMatrix<double, Eigen::RowMajor> braket_remove(const SparseMatrix<double, Eigen::RowMajor> &X, const MatrixXi &ref);

	SparseMatrix<double> RowBaseToColBase(SparseMatrix<double, Eigen::RowMajor> X);

	SparseMatrix<double, Eigen::RowMajor> ColBaseToRowBase(SparseMatrix<double> X);

	Mesh StreamToMesh(double vertices[], int VerticesLen, int triangles, int trianglesLen);

	class KnnSearch
	{
	public:
		KnnSearch(MatrixXd vertices);
		MatrixXi Kneighbors(MatrixXd vertices, int K =1,double threhold = __THREHOLD);// get K neighbors
	private:
		KDTree kdtree;
	};
	
	mesh_v_and_f get__head(::std::vector<Eigen::Vector3f> vertice, ::std::vector<::std::array<int, 3>> triangles);

	Mesh Fit(Mesh x);
	Mesh Fit(Mesh face, Mesh head);


}

extern "C" {

	STITCH_API void __C__mesh__test();
	STITCH_API void __C__mesh__test2();
	STITCH_API void neck_test();


}

