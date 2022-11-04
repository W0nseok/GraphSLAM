#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>

class PoseEdge
{
	public:
		PoseEdge(int id_from_, int id_to_, Eigen::MatrixXd mean_, Eigen::MatrixXd infm_);

		int id_from;
		int id_to;
		Eigen::MatrixXd mean;
		Eigen::MatrixXd infm;
};

class Posegraph
{
	public:
		Posegraph();

		Eigen::MatrixXd nodes;
		std::vector<PoseEdge> edges;
		Eigen::MatrixXd landmark_nodes;
		std::vector<PoseEdge> landmark_edges;
		Eigen::MatrixXd H;
		Eigen::MatrixXd b;
		// edge를 pose edge로 받아오는 코드

		void readGraph(std::string vfile, std::string efile);

		void optimize();

		void linearize();

		void solve();

		void print();

		int num_of_nodes;
		int num_of_l_nodes;
};