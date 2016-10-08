#include "brg.h"

namespace brg
{
	BRG_API void generate_graph(int N, int M, unsigned seed, Edge** result)
	{
		const int kN = N;
		const int kM = M;

		std::default_random_engine gen(seed);
		std::uniform_real_distribution<float> distr(0.0, 1.0);

		std::vector<std::list<int> > adj_list_nm_1(kN * kM);
		std::vector<int> ranks_nm_1(kN * kM);
		for (int i = 1; i <= kN * kM; ++i) {
			ranks_nm_1[i - 1]++;

			float total_rank = 2 * i - 1;
			float prob = distr(gen);

			float acc = 0;
			bool edge = false;
			for (int j = 0; j < i; ++j) {
				if (acc <= prob && prob < acc + ranks_nm_1[j] / total_rank) {
					edge = true;
					ranks_nm_1[j]++;
					adj_list_nm_1[j].push_back(i - 1);
					if (i - 1 != j) adj_list_nm_1[i - 1].push_back(j);
					break;
				}
				acc += ranks_nm_1[j] / total_rank;
			}
			if (!edge) {
				ranks_nm_1[i - 1]++;
				adj_list_nm_1[i - 1].push_back(i - 1);
			}
		}

		std::vector<std::list<int> > adj_list_n_m(kN);
		for (int i = 0; i < kN; ++i) {
			for (int j = i * kM; j < (i + 1) * kM; ++j) {
				for (std::list<int>::const_iterator it = adj_list_nm_1[j].begin(); it != adj_list_nm_1[j].end(); ++it) {
					if (*it >= j) adj_list_n_m[i].push_back(*it / kM);
					if (*it >= (i + 1) * kM) adj_list_n_m[*it / kM].push_back(i);
				}
			}
		}

		Edge* edges = new Edge[kN * kM];
		*result = edges;
		for (int i = 0, j = 0; i < kN; ++i) {
			for (std::list<int>::const_iterator it = adj_list_n_m[i].begin(); it != adj_list_n_m[i].end(); ++it) {
				if (*it >= i) edges[j++] = {i, *it};
			}
		}
	}
}