#include "posegraph.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

Posegraph::Posegraph()
{
	num_of_nodes = 0;
	num_of_l_nodes = 0;
}

PoseEdge::PoseEdge(int id_from_, int id_to_, Eigen::MatrixXd mean_, Eigen::MatrixXd infm_)
{
	id_from = id_from_;
	id_to = id_to_;
	mean = mean_;
	infm = infm_;
}

std::vector<std::string> split(std::string input, char delimiter)
{
	std::vector<std::string> answer;
	std::stringstream ss(input);
	std::string temp;

	while (getline(ss, temp, delimiter))
	{
		answer.push_back(temp);
	}

	return answer;
}

void Posegraph::readGraph(std::string input_vfile, std::string input_efile)
{
	std::ifstream vfile(input_vfile);
	std::ifstream efile(input_efile);
	Eigen::Matrix3d infm;
	Eigen::Matrix2d linfm;
	Eigen::MatrixXd mean(3, 1);
	Eigen::MatrixXd lmean(2, 1);
	int id_from, id_to;
	std::string line;

	std::vector<Eigen::MatrixXd> temp_node;
	std::vector<Eigen::MatrixXd> temp_lnode;

	if (vfile.is_open())
	{
		while (vfile)
		{
			getline(vfile, line);
			std::cout << line << std::endl;
			if (!line.empty())
			{
				std::vector<std::string> buffer = split(line, ' ');

				if (buffer[0] == "V2")
				{

					Eigen::MatrixXd node(1, 3);

					node(0, 0) = stod(buffer[2]);
					node(0, 1) = stod(buffer[3]);
					node(0, 2) = stod(buffer[4]);

					temp_node.push_back(node);

					num_of_nodes++;
				}
				else if (buffer[0] == "L")
				{
					Eigen::MatrixXd lnode(1, 2);

					lnode(0, 0) = stod(buffer[2]);
					lnode(0, 1) = stod(buffer[3]);

					temp_lnode.push_back(lnode);

					num_of_l_nodes++;
				}
			}
		}
	}

	nodes = Eigen::MatrixXd::Zero(num_of_nodes, 3);

	for (int i = 0; i < num_of_nodes; i++)
	{
		nodes.block(i, 0, 1, 3) = temp_node[i];
	}

	landmark_nodes = Eigen::MatrixXd::Zero(num_of_l_nodes, 2);

	for (int i = 0; i < num_of_l_nodes; i++)
	{
		landmark_nodes.block(i, 0, 1, 2) = temp_lnode[i];
	}

	if (efile.is_open())
	{
		while (efile)
		{
			getline(efile, line);
			if (!line.empty())
			{
				std::vector<std::string> buffer = split(line, ' ');

				if (buffer[0] == "EDGE2")
				{

					mean(0, 0) = stod(buffer[3]);
					mean(1, 0) = stod(buffer[4]);
					mean(2, 0) = stod(buffer[5]);

					id_from = stoi(buffer[1]);
					id_to = stoi(buffer[2]);

					infm(0, 0) = stod(buffer[6]);
					infm(1, 0) = infm(0, 1) = stod(buffer[7]);
					infm(1, 1) = stod(buffer[8]);
					infm(2, 2) = stod(buffer[9]);
					infm(0, 2) = infm(0, 2) = stod(buffer[10]);
					infm(1, 2) = infm(2, 1) = stod(buffer[11]);

					//				PoseEdge edge = PoseEdge(id_from, id_to, mean, infm);
					PoseEdge edge = PoseEdge(id_to, id_from, mean, infm);

					edges.push_back(edge);
				}

				else if (buffer[0] == "LE")
				{
					lmean(0, 0) = stod(buffer[3]);
					lmean(1, 0) = stod(buffer[4]);

					id_from = stoi(buffer[1]); // robot pose
					id_to = stoi(buffer[2]); // landmark number

					linfm(0, 0) = stod(buffer[5]);
					linfm(1, 0) = linfm(0, 1) = stod(buffer[6]);
					linfm(1, 1) = stod(buffer[7]);

					//				PoseEdge edge = PoseEdge(id_from, id_to, mean, infm);
					PoseEdge ledge = PoseEdge(id_from, id_to, lmean, linfm);

					landmark_edges.push_back(ledge);
				}
			}
		}
	}
}


void Posegraph::optimize()
{
	int iter = 5;
	for (int i = 0; i < iter; i++)
	{
		std::cout << "optimze start" << std::endl;

		H = Eigen::MatrixXd::Zero(3 * num_of_nodes + 2 * num_of_l_nodes, 3 * num_of_nodes + 2 * num_of_l_nodes);
		b = Eigen::MatrixXd::Zero(3 * num_of_nodes + 2 * num_of_l_nodes, 1);

		Posegraph::linearize();

		Posegraph::solve();

		std::cout << "print result" << std::endl;

		Posegraph::print();
	}
}

void Posegraph::linearize()
{
	std::cout << "linearize start" << std::endl;

	for (int i = 0; i < edges.size(); i++)
	{
		PoseEdge edge = edges[i];

		int i_node = edge.id_from;
		int j_node = edge.id_to;

		Eigen::MatrixXd T_z(3, 3);
		T_z = v2t(edge.mean);
		Eigen::MatrixXd Omega(3, 3);
		Omega = edge.infm;

		Eigen::MatrixXd v_i(3, 1);
		Eigen::MatrixXd v_j(3, 1);

		v_i = nodes.row(i_node).transpose(); // 1 x 3이어서 transpose
		v_j = nodes.row(j_node).transpose();

		Eigen::MatrixXd T_i(3, 3);
		Eigen::MatrixXd T_j(3, 3);

		T_i = v2t(v_i);
		T_j = v2t(v_j);

		Eigen::MatrixXd R_i(2, 2);
		Eigen::MatrixXd R_z(2, 2);

		R_i = T_i.block(0, 0, 2, 2);
		R_z = T_z.block(0, 0, 2, 2);

		double si = sin(v_i(2, 0));
		double co = cos(v_i(2, 0));

		Eigen::MatrixXd dR_i(2, 2);
		dR_i << -si, -co, co, -si;

		Eigen::MatrixXd dt_ij(2, 1);
		dt_ij << v_j(0, 0) - v_i(0, 0), v_j(1, 0) - v_i(1, 0);

		Eigen::MatrixXd A(3, 3);
		A.block(0, 0, 2, 2) = -R_z.transpose() * R_i.transpose();
		A.block(0, 2, 2, 1) = R_z.transpose() * dR_i.transpose() * dt_ij;
		A(2, 0) = 0;
		A(2, 1) = 0;
		A(2, 2) = -1;

		Eigen::MatrixXd B(3, 3);
		B.block(0, 0, 2, 2) = R_z.transpose() * R_i.transpose();
		B(0, 2) = 0;
		B(1, 2) = 0;
		B(2, 0) = 0;
		B(2, 1) = 0;
		B(2, 2) = 1;

		Eigen::MatrixXd e(3, 1);
		e = t2v(T_z.inverse() * T_i.inverse() * T_j);

		Eigen::MatrixXd H_ii(3, 3);
		Eigen::MatrixXd H_ij(3, 3);
		Eigen::MatrixXd H_jj(3, 3);
		H_ii = A.transpose() * Omega * A;
		H_ij = A.transpose() * Omega * B;
		H_jj = B.transpose() * Omega * B;

		Eigen::MatrixXd b_i(3, 1);
		Eigen::MatrixXd b_j(3, 1);
		b_i = -A.transpose() * Omega * e;
		b_j = -B.transpose() * Omega * e;

		H.block(3 * i_node, 3 * i_node, 3, 3) += H_ii;
		H.block(3 * i_node, 3 * j_node, 3, 3) += H_ij;
		H.block(3 * j_node, 3 * i_node, 3, 3) += H_ij.transpose();
		H.block(3 * j_node, 3 * j_node, 3, 3) += H_jj;

		b.block(3 * i_node, 0, 3, 1) += b_i;
		b.block(3 * j_node, 0, 3, 1) += b_j;
	}

	for (int i = 0; i < landmark_edges.size(); i++)
	{
		PoseEdge landmark_edge = landmark_edges[i];

		int i_node = landmark_edge.id_from;
		int j_node = landmark_edge.id_to;

		Eigen::MatrixXd T_z(2, 1);
		T_z(0, 0) = landmark_edge.mean(0, 0) * cos(landmark_edge.mean(1, 0));
		T_z(1, 0) = landmark_edge.mean(0, 0) * sin(landmark_edge.mean(1, 0));
		Eigen::MatrixXd Omega(2, 2);
		Omega = landmark_edge.infm;
		
		Eigen::MatrixXd v_i(3, 1); // robot (x, y, theta)
		Eigen::MatrixXd v_j(2, 1); // landmark (x, y)

		v_i = nodes.row(i_node).transpose(); // 1 x 3이어서 transpose
		v_j = landmark_nodes.row(j_node).transpose();

		Eigen::MatrixXd T_i(3, 3);
		// Eigen::MatrixXd T_j(3, 3);

		T_i = v2t(v_i);
		// T_j = v2t(v_j);

		Eigen::MatrixXd R_i(2, 2);
		// Eigen::MatrixXd R_z(2, 2);

		R_i = T_i.block(0, 0, 2, 2);
		// R_z = T_z.block(0, 0, 2, 2);

		double si = sin(v_i(2, 0));
		double co = cos(v_i(2, 0));

		Eigen::MatrixXd dR_i(2, 2);
		dR_i << -si, -co, co, -si;

		Eigen::MatrixXd dt_ij(2, 1);
		dt_ij << v_j(0, 0) - v_i(0, 0), v_j(1, 0) - v_i(1, 0);

		Eigen::MatrixXd A(2, 3);
		A.block(0, 0, 2, 2) = -R_i.transpose();
		A.block(0, 2, 2, 1) = dR_i.transpose() * dt_ij;
		
		Eigen::MatrixXd B(2, 2);
		B.block(0, 0, 2, 2) = R_i.transpose();
		
		Eigen::MatrixXd e(2, 1);
		e = R_i.transpose() * dt_ij - T_z;

		Eigen::MatrixXd H_ii(3, 3);
		Eigen::MatrixXd H_ij(3, 2);
		Eigen::MatrixXd H_jj(2, 2);
		H_ii = A.transpose() * Omega * A;
		H_ij = A.transpose() * Omega * B;
		H_jj = B.transpose() * Omega * B;

		Eigen::MatrixXd b_i(3, 1);
		Eigen::MatrixXd b_j(2, 1);
		b_i = -A.transpose() * Omega * e;
		b_j = -B.transpose() * Omega * e;

		H.block(3 * i_node, 3 * i_node, 3, 3) += H_ii;
		H.block(3 * i_node, 3 * num_of_nodes + 2 * j_node, 3, 2) += H_ij;
		H.block(3 * num_of_nodes + 2 * j_node, 3 * i_node, 2, 3) += H_ij.transpose();
		H.block(3 * num_of_nodes + 2 * j_node, 3 * num_of_nodes + 2 * j_node, 2, 2) += H_jj;

		b.block(3 * i_node, 0, 3, 1) += b_i;
		b.block(3 * num_of_nodes + 2 * j_node, 0, 2, 1) += b_j;
	}
}

void Posegraph::solve()
{
	H.block(0, 0, 3, 3) += Eigen::Matrix3d::Identity();

	Eigen::SparseMatrix<double> H_sparse = H.sparseView();
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(H_sparse);
	Eigen::MatrixXd temp_dx(3 * num_of_nodes + 2 * num_of_l_nodes, 1);
	Eigen::MatrixXd dx(num_of_nodes, 3);
	
	temp_dx = chol.solve(b);
	temp_dx(0, 0) = 0;
	temp_dx(1, 0) = 0;
	temp_dx(2, 0) = 0;
	//dx.block(0, 0, 3, 1) << 0, 0, 0;
	

	if (num_of_l_nodes == 0)
	{
		temp_dx = temp_dx.reshaped(3, num_of_nodes);

		dx = temp_dx.transpose();

		std::cout << "before update: \n" << nodes << std::endl;

		std::cout << "dx: \n" << dx << std::endl;

		nodes += dx;

		std::cout << "after update: \n" << nodes << std::endl;
	}
	else
	{

		Eigen::MatrixXd l_dx(num_of_l_nodes, 2);
		Eigen::MatrixXd temp_l_dx(2 * num_of_l_nodes, 1);
		temp_l_dx = temp_dx.block(3 * num_of_nodes, 0, 2 * num_of_l_nodes, 1);

		temp_l_dx = temp_l_dx.reshaped(2, num_of_l_nodes);
		l_dx = temp_l_dx.transpose();

		temp_dx = temp_dx.block(0, 0, 3 * num_of_nodes, 1);

		temp_dx = temp_dx.reshaped(3, num_of_nodes);

		dx = temp_dx.transpose();


		std::cout << "before update: \n" << nodes << std::endl;
		std::cout << "before L-update: \n" << landmark_nodes << std::endl;
		std::cout << "dx: \n" << dx << std::endl;
		std::cout << "landmark_dx: \n" << l_dx << std::endl;

		nodes += dx;

		landmark_nodes += l_dx;

		std::cout << "after update: \n" << nodes << std::endl;
		std::cout << "after L-update: \n" << landmark_nodes << std::endl;

	}
}

void Posegraph::print()
{
	int print_num = num_of_nodes;
	int print_l_num = num_of_l_nodes;
	//int print_num = 100;
	std::cout << "pose" << std::endl;
	for (int i = 0; i < print_num; i++)
	{
		std::cout << nodes.row(i) << std::endl;
	}
	if (num_of_l_nodes != 0)
	{
		std::cout << "landmark" << std::endl;
		for (int i = 0; i < print_l_num; i++)
		{
			std::cout << landmark_nodes.row(i) << std::endl;
		}
	}

	bool print = false;

	if (print)
	{
		std::ofstream writeFile;
		writeFile.open("pose_0.01.txt");
		Eigen::IOFormat EigenFormat(6);
		for (int i = 0; i < print_num; i++)
		{
			writeFile << nodes.row(i).format(EigenFormat);
			writeFile << "\n";
		}
		writeFile.close();

		writeFile.open("landmark_0.01.txt");
		for (int i = 0; i < print_l_num; i++)
		{
			writeFile << landmark_nodes.row(i).format(EigenFormat);
			writeFile << "\n";
		}
		writeFile.close();
	}
	
}

