#include "ParallelLeiden.hpp"
#include <networkit/io/METISGraphReader.hpp>
#include <fstream>
#include <networkit/io/SNAPGraphReader.hpp>


int main() {
    Graph g = SNAPGraphReader().read("../LFR0,7.el");
    ParallelLeiden pl(g);
    pl.setMaxThreads(8);
    auto tmr = Aux::Timer();
    tmr.start();
    pl.run();
    cout << tmr.elapsedMilliseconds() / 1000.0 << " " << Modularity().getQuality(pl.getPartition(), g) << endl;
    return 0;
}

