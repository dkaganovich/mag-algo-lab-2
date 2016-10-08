#pragma once

#include <vector>
#include <list>
#include <random>

#ifdef BRG_EXPORTS
#define BRG_API __declspec(dllexport) 
#else
#define BRG_API __declspec(dllimport) 
#endif

namespace brg 
{

	struct Edge {
		int v0;
		int v1;
	};

	BRG_API void generate_graph(int N, int M, unsigned seed, Edge** result);

}