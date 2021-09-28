#ifndef ParallelLeidenHPP
#define ParallelLeidenHPP

//#define PLDEBUG // Remove for release
//#define PLFINE
//#define PLBENCHMARK
#include <omp.h>

#pragma GCC system_header // Ain't nobody got time for Library GCC warnings. Shut up.

#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <iostream>
#include <networkit/auxiliary/Parallel.hpp>
#include <mutex>
#include <networkit/auxiliary/Timer.hpp>
#include <atomic>
#include <condition_variable>
#include <thread>
#include <networkit/auxiliary/Parallelism.hpp>
#include <iomanip>
#include <chrono>
#include <networkit/community/Modularity.hpp>
#include <filesystem>
#include <fstream>
#include <utility>
#include <iostream>


#ifdef PLDEBUG

#define PLPRINT(x) cout << x << "\n";
#define PLIF(x) x;
#else
#define PLPRINT(x)
#define PLIF(x)

#endif

#ifdef PLFINE
#define PLIFF(x) x;
#define PLTRACE(x) cout << x << "\n";
#else
#define PLTRACE(x)
#define PLIFF(x)
#endif

using namespace NetworKit;
using namespace std;

namespace NetworKit {

    class ParallelLeiden final : public CommunityDetectionAlgorithm {
    public:
        explicit ParallelLeiden(const Graph &graph, int iterations = 3, bool randomize = true, double gamma = 1);

        ParallelLeiden(const Graph &graph, Partition basePartition, int iterations = 3, bool randomize = true, double gamma = 1);

        void run() override;

        enum parallelModes {
            SEQUENTIAL, PARALLEL, EXPERIMENTAL
        };

        void setParallelism(parallelModes localMoving, parallelModes refinement) {
            this->plMove = localMoving;
            this->plRefine = refinement;
        }

        void setMaxThreads(int count) {
            this->tC = count;
        }

        void setGraphName(string name) {
            this->graphName = std::move(name);
        }


    private :

        [[nodiscard]] inline double modularityDelta(double cutD, double degreeV, double volD) const {
            return cutD - gamma * degreeV * volD * inverseGraphVolume;
        };

        [[nodiscard]] inline double modularityThreshold(double cutC, double volC, double degreeV) const {
            return cutC - gamma * (volC - degreeV) * degreeV * inverseGraphVolume;
        }

        static inline void lockLowerFirst(index a, index b, vector<mutex> &locks){
            if(a < b){
                locks[a].lock();
                locks[b].lock();
            }
            else{
                locks[b].lock();
                locks[a].lock();
            }
        }

        void flattenPartition();

        void calculateVolumes(const Graph &graph);


        //Sequential
        void moveNodesFast(const Graph &graph);

        Partition refineAndMerge(const Graph &graph);

        //Parallel
        void parallelMove(const Graph &graph);

        Partition parallelRefine(const Graph &graph);

        //Experimental
        void experimentalMove(const Graph &graph);

        Partition experimentalRefine(const Graph &graph);

        static string enumToString(parallelModes mode) {
            switch (mode) {
                case SEQUENTIAL:
                    return "SQ";
                case PARALLEL:
                    return "PL";
                case EXPERIMENTAL:
                    return "EXP";
            }
        }

        double inverseGraphVolume; // 1/vol(V)

        vector<double> communityVolumes;

        vector<vector<node>> mappings;

        double gamma;

        bool changed = true;

        int numberOfIterations;

        Aux::SignalHandler handler;

        int tC = omp_get_max_threads();

        string graphName;

        int rounds = 0;

        parallelModes plMove = PARALLEL;

        parallelModes plRefine = PARALLEL;

        bool random;
    };

}

#endif //ParallelLeidenHPP
