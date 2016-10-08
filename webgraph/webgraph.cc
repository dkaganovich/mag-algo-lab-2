#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <chrono>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/format.hpp>


#include <make_graph.h>
#include <brg.h>

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>

namespace {
	// convenience types
	typedef int64_t VertexId;

	namespace bst = boost;

	typedef bst::property<bst::vertex_rank_t, int, bst::property<bst::vertex_name_t, int64_t> > VertexProps;
	typedef bst::property<bst::edge_weight_t, int> EdgeProps;
	typedef bst::adjacency_list<bst::vecS, bst::vecS, bst::undirectedS, VertexProps, EdgeProps> Graph;

	typedef bst::property_map<Graph, bst::vertex_rank_t>::type RankProps;
	typedef bst::property_map<Graph, bst::edge_weight_t>::type WeightProps;

	typedef bst::graph_traits<Graph>::vertex_descriptor VertexDesc;
	typedef bst::graph_traits<Graph>::edge_descriptor EdgeDesc;

	typedef std::unordered_map<VertexId, VertexDesc> VertexStore;

	typedef std::unordered_map<int, float> RankHist;

	typedef bst::exterior_vertex_property<Graph, int> DistanceProps;
	typedef DistanceProps::matrix_type DistanceMx;
	typedef DistanceProps::matrix_map_type DistanceMxMap;

	void generate_graph(const std::string& algo, int N, int M, Graph& g, RankProps& rank_props, WeightProps& weight_props) 
	{
		VertexStore verts_store;
		for (int i = 0; i < N; ++i) {
			VertexDesc v = add_vertex(g);
			verts_store[i] = v;
			rank_props[v] = 0;
		}

		if (algo == "--skg") {
			if ((N & (N - 1)) != 0) {
				std::cerr << "For --skg mode N must be a power of 2" << std::endl;
				exit(1);
			}

			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

			int64_t nedges;
			packed_edge* result;
			make_graph((int)log2f(N), (int64_t)N * M, (int64_t)seed, (int64_t)seed + 1, &nedges, &result);

			for (int i = 0; i < nedges; ++i) {
				VertexDesc v0 = verts_store[get_v0_from_edge(&result[i])];
				VertexDesc v1 = verts_store[get_v1_from_edge(&result[i])];
				std::pair<EdgeDesc, bool> e = add_edge(verts_store[v0], verts_store[v1], g);
				weight_props[e.first] = 1;
				rank_props[v0]++;
				rank_props[v1]++;
			}

			free(result);
		}
		else {// --brg
			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

			brg::Edge* result;
			brg::generate_graph(N, M, seed, &result);

			for (int i = 0; i < N * M; ++i) {
				VertexDesc v0 = verts_store[result[i].v0];
				VertexDesc v1 = verts_store[result[i].v1];
				std::pair<EdgeDesc, bool> e = add_edge(verts_store[v0], verts_store[v1], g);
				weight_props[e.first] = 1;
				rank_props[v0]++;
				rank_props[v1]++;
			}

			delete[] result;
		}
	}
}

int main(int argc, char* argv[])
{
	std::string algo;
	if (argc >= 2) algo = std::string(argv[1]);
	if (algo != "--skg" && algo != "--brg") {
		std::cerr << ("Unknown algorithm") << std::endl;
		exit(1);
	}
	int N = 4096;// num vertices
	if (argc >= 3) N = atoi(argv[2]);
	int M = 16;// edge factor
	if (argc >= 4) M = atoi(argv[3]);
	int numTrials = 16;
	if (argc >= 5) numTrials = atoi(argv[4]);

	const std::string kSmallWorldFname = bst::str(bst::format("%1%_small_world.csv") % algo);
	const std::string kNumCompsFname = bst::str(bst::format("%1%_num_comps.csv") % algo);
	const std::string kRankHistFname = bst::str(bst::format("%1%_rank_hist.csv") % algo);

	std::ofstream ofs;

	RankHist rank_hist;
	for (int i = 0; i < numTrials; ++i)
	{
		std::cout << "Trial " << i << ": started: " << 1. * std::clock() / CLOCKS_PER_SEC << std::endl;

		Graph g;
		RankProps rank_props = get(bst::vertex_rank, g);
		WeightProps weight_props = get(bst::edge_weight, g);

		generate_graph(algo, N, M, g, rank_props, weight_props);

		// rank stats
		bst::graph_traits<Graph>::vertex_iterator vi, vi_end;
		for (bst::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
			rank_hist[rank_props[*vi]] += 1. / N;
		}

		// small world test
		DistanceMx dists(N);
		DistanceMxMap dists_map(dists, g);
		johnson_all_pairs_shortest_paths(g, dists_map, bst::weight_map(weight_props));// no negative cycles

		float diam = 0;
		float dist = 0;
		for (std::size_t i = 0; i < N; ++i) {
			float dist_ = 0;
			int cnt = 0;
			for (std::size_t j = 0; j < N; ++j) {
				if (dists[i][j] != 0 && dists[i][j] != std::numeric_limits<int>::max()) {
					if (diam < dists[i][j]) diam = dists[i][j];
					dist_ += dists[i][j];
					cnt++;
				}
			}
			if (cnt > 0)
				dist += dist_ / cnt;
		}
		// output
		ofs.open(kSmallWorldFname, std::ofstream::out | std::ofstream::app);
		ofs << diam << "," << dist / N << std::endl;
		ofs.close();

		// connected components test
		std::vector<int> component(N);
		int num_comps = connected_components(g, &component[0]);
		// output
		ofs.open(kNumCompsFname, std::ofstream::out | std::ofstream::app);
		ofs << num_comps << std::endl;
		ofs.close();

		std::cout << "Trial " << i << ": finished: " << 1. * std::clock() / CLOCKS_PER_SEC << std::endl;
	}
	// output rank stats
	ofs.open(kRankHistFname, std::ofstream::out);
	for (auto it = rank_hist.begin(); it != rank_hist.end(); ++it) {
		ofs << it->first << "," << 1. * it->second / numTrials << std::endl;
	}
	ofs.close();
}