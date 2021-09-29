#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"

#include "ParallelLeiden.hpp"

#define WORKING_SIZE 1000

namespace NetworKit {

    /// Leiden Algorithm
    /// \param graph The graph to run the leiden algorithm on
    /// \param iterations How many times the leiden algorithl should be run. Default is 3, stops early if the partition doesn't change anymore
    /// \param gamma Resolution parameter between 0(one community) and 2m(singletons). Default 1
    ParallelLeiden::ParallelLeiden(const Graph &graph, int iterations, bool randomize, double gamma) :
            CommunityDetectionAlgorithm(graph), gamma(gamma), numberOfIterations(iterations), random(randomize) {
        this->result = Partition(graph.numberOfNodes());
        this->result.allToSingletons();
    }

    ParallelLeiden::ParallelLeiden(const Graph &graph, Partition basePartition, int iterations, bool randomize, double gamma) :
            CommunityDetectionAlgorithm(graph, std::move(basePartition)), gamma(gamma), numberOfIterations(iterations), random(randomize) {}

    void ParallelLeiden::run() {
        if (plMove == EXPERIMENTAL || plRefine == EXPERIMENTAL)
            cerr << "You're using an experimental function that has not been tested thoroughly. Unless you know what you're doing you should be using the standard functions." << endl;

        Aux::setNumberOfThreads(tC);
        cout << "Running leiden Algorithm" << endl;
        double mv = 0.0;
        double ref = 0.0;
        double agg = 0.0;
        auto totalTime = Aux::Timer();
        totalTime.start();
        do {            // Leiden recursion
            PLPRINT("Leiden rekursion " << numberOfIterations << " left")
            numberOfIterations--;
            changed = false;
            bool done;                           //have we moved any nodes? If yes, start local moving phase again. Else stop
            const Graph *currentGraph = G;
            Graph coarse;
            Partition refined;
            calculateVolumes(*currentGraph);
            do {
                auto tM = Aux::Timer();
                tM.start();
                switch (plMove) {
                    case SEQUENTIAL:
                        moveNodesFast(*currentGraph);
                        break;
                    case PARALLEL:
                        parallelMove(*currentGraph);
                        break;
                    case EXPERIMENTAL:
                        experimentalMove(*currentGraph);
                        break;
                    default:
                        throw;
                }
                mv += tM.elapsedMilliseconds() / 1000.0;
                done = currentGraph->numberOfNodes() == result.numberOfSubsets();            // If each community consists of exactly one node we're done, i.e. when |V(G)| = |P|
                if (!done) {
                    auto tR = Aux::Timer();
                    tR.start();
                    switch (plRefine) {
                        case SEQUENTIAL:
                            refined = refineAndMerge(*currentGraph);
                            break;
                        case PARALLEL:
                            refined = parallelRefine(*currentGraph);
                            break;
                        case EXPERIMENTAL:
                            refined = experimentalRefine(*currentGraph);
                            break;
                        default:
                            throw;
                    }

                    ref += tR.elapsedMilliseconds() / 1000.0;
                    auto tPPC = Aux::Timer();
                    tPPC.start();

                    ParallelPartitionCoarsening ppc(*currentGraph, refined);                 //Aggregate graph
                    ppc.run();
                    auto temp = move(ppc.getCoarseGraph());
                    agg += tPPC.elapsedMilliseconds() / 1000.0;

                    auto map = std::move(ppc.getFineToCoarseNodeMapping());

                    Partition p(temp.numberOfNodes());                                    // "Maintain Partition" : add every coarse Node to the community its fine Nodes were in
                    p.setUpperBound(
                            result.upperBound());                             // this isn't needed for louvain but with leiden 2 coarse Nodes can belong to the same community

                    currentGraph->parallelForNodes([&](node Node) {
                        p[map[Node]] = result[Node];
                    });

                    mappings.emplace_back(std::move(map));
                    result = std::move(p);
                    swap(temp, coarse);
                    currentGraph = &coarse;
                }
                PLPRINT("---------------------------------------------------")
            } while (!done);
            flattenPartition();
            PLPRINT("Leiden done. Modularity: " << setprecision(8) << Modularity().getQuality(result, *G))
        } while (changed && numberOfIterations > 0);
        hasRun = true;
    }

    void ParallelLeiden::calculateVolumes(const Graph &graph) {
#ifdef PLFINE
        auto timer = Aux::Timer();
        timer.start();
#endif
        // thread safe reduction. Avoid atomic calculation of total graph volume for unweighted graphs. Vol(G) is then 2*|E|
        communityVolumes.clear();
        communityVolumes.resize(result.upperBound() + 10000);
        if (graph.isWeighted()) {
            vector<double> threadVolumes(omp_get_max_threads());
            graph.parallelForNodes([&](node a) {
                {
                    edgeweight ew = graph.weightedDegree(a, true);
#pragma omp atomic
                    communityVolumes[result[a]] += ew;
                    threadVolumes[omp_get_thread_num()] += ew;
                }
            });
            for (const auto vol: threadVolumes) {
                inverseGraphVolume += vol;
            }
            inverseGraphVolume = 1 / inverseGraphVolume;
        } else {
            inverseGraphVolume = 1.0 / (2 * graph.numberOfEdges());
            graph.parallelForNodes([&](node a) {
                {
#pragma omp atomic
                    communityVolumes[result[a]] += graph.weightedDegree(a, true);
                }
            });
        }
        PLTRACE("Calculating Volumes took " << timer.elapsedMilliseconds() / 1000.0 << "s")
    }

    void ParallelLeiden::flattenPartition() {
#ifdef PLFINE
        auto timer = Aux::Timer();
        timer.start();
#endif
        if (mappings.empty()) {
            return;
        }
        Partition flattenedPartition(G->numberOfNodes());           // Create a new partition with size |V(G)| (the fine/bigger Graph)
        flattenedPartition.setUpperBound(result.upperBound());
        int i = mappings.size() - 1;
        vector<node> &lower = mappings[i--];
        while (i >= 0) {                                   // iteratively "resolve" (i.e compose) mappings. Let "lower" be a mapping thats below "higher" in the hierarchy (i.e. of a later aggregation)
            vector<node> &upper = mappings[i--];                     // If higher[index] = z and lower[z] = x then: higher[index] = x
            for (auto &idx: upper) {
                idx = lower[idx];
            }
            lower = upper;
        }
        G->parallelForNodes([&](node a) {
            flattenedPartition[a] = lower[a];
        });
        flattenedPartition.compact(true);
        result = flattenedPartition;
        mappings.clear();
        PLTRACE("Flattening partition took " << timer.elapsedMilliseconds() / 1000.0 << "s")
    }

    void ParallelLeiden::parallelMove(const Graph &graph) {
#ifdef PLFINE
        cout << "Local Moving : " << graph.numberOfNodes() << " Nodes " << endl;
        vector<long> moved(omp_get_max_threads(), 0);
        vector<long> totalNodesPerThread(omp_get_max_threads(), 0);
        atomic_int singleton = 0;
#endif
        vector<atomic_bool> inQueue(graph.upperNodeIdBound());      // Only insert nodes to the queue when they're not already in it.
        queue<vector<node>> queue;
        mutex qlock;    // queue lock
        condition_variable workAvailable; // waiting/notifying for new Nodes

        atomic_bool resize = false;
        atomic_int waitingForResize = 0;
        atomic_int waitingForNodes = 0;

        vector<int> order;
        int tshare;
        int tcount;
        uint64_t vectorSize = communityVolumes.capacity();
        atomic_int upperBound = result.upperBound();
#pragma omp parallel
        {
#pragma omp single
            {
                tcount = omp_get_num_threads();
                order.resize(tcount);
                for (int i = 0; i < tcount; i++) {
                    order[i] = i;
                }
                if (random)
                    shuffle(order.begin(), order.end(), Aux::Random::getURNG());
                tshare = 1 + graph.upperNodeIdBound() / tcount;
            }
            auto &mt = Aux::Random::getURNG();
            vector<node> currentNodes;
            currentNodes.reserve(tshare);
            vector<node> newNodes;
            newNodes.reserve(WORKING_SIZE);
            vector<double> cutWeights(communityVolumes.capacity());      // cutWeight[Community] returns cut of Node to Community
            vector<index> pointers;
            int start = tshare * order[omp_get_thread_num()];
            int end = (1 + order[omp_get_thread_num()]) * tshare;

            for (int i = start; i < end; i++) {
                if (graph.hasNode(i)) {
                    currentNodes.push_back(i);
                    atomic_init(&inQueue[i], true);
                }
            }
            if (random)
                shuffle(currentNodes.begin(), currentNodes.end(), mt);
#pragma omp barrier
            do {
                for (node Node: currentNodes) {
                    if (resize) {       // This will probably never happen. community IDs may be larger than initial size of vector if new IDs are assigned to singletons
                        waitingForResize++;
                        while (resize) {
                            this_thread::yield();
                        }
                        waitingForResize--;
                    }
                    cutWeights.resize(vectorSize);
                    assert(inQueue[Node]);
                    index currentCommunity = result[Node];                    // The node's current communityID
                    double maxDelta = numeric_limits<double>::lowest();       // remember best delta so far ..
                    index bestCommunity = none;                               // and the corresponding communityID
                    double degree = 0;                                      // weighted degree of the Node
                    for (auto z: pointers) {                               // Reset the clearlist : Set all cutweights to 0 and clear the pointer vector
                        cutWeights[z] = 0;
                    }
                    pointers.clear();

                    graph.forNeighborsOf(Node, [&](node neighbor, edgeweight ew) {
                        index neighborCommunity = result[neighbor];
                        if (cutWeights[neighborCommunity] == 0) {
                            pointers.push_back(neighborCommunity);
                        }
                        if (Node == neighbor) {
                            degree += ew;
                        } else {
                            cutWeights[neighborCommunity] += ew;
                        }
                        degree += ew;                                       // keep track of the nodes degree, since we're already iterating over its neighbors. Loops count twice
                    });

                    if (pointers.empty())
                        continue;

                    for (auto community: pointers) {                                                                 // For all neighbor communities, determine the modularity Delta for moving the node to the community
                        if (community != currentCommunity) {                                                            // "Moving" a node to its current community is pointless
                            double delta;
                            delta = modularityDelta(cutWeights[community], degree, communityVolumes[community]);
                            if (delta > maxDelta) {                                                                 // Keep track of the best delta and the corresponding communityID
                                maxDelta = delta;
                                bestCommunity = community;
                            }
                        }
                    }
                    double modThreshold = modularityThreshold(cutWeights[currentCommunity], communityVolumes[currentCommunity], degree);

                    if (0 > modThreshold || maxDelta > modThreshold) {
                        PLIFF(moved[omp_get_thread_num()]++)
                        if (0 > maxDelta) {                         // move node to empty community
                            PLIFF(singleton++)
                            bestCommunity = upperBound++;
                            if (bestCommunity >= communityVolumes.capacity()) {                 // Chances are this will never happen. Ever.
                                bool expected = false;
                                if (resize.compare_exchange_strong(expected, true)) {
                                    vectorSize += 10000;
                                    cutWeights.resize(vectorSize);
                                    while (waitingForResize < tcount - 1) {
                                        this_thread::yield();
                                    }
                                    communityVolumes.resize(vectorSize);
                                    expected = true;
                                    resize.compare_exchange_strong(expected, false);
                                } else {
                                    waitingForResize++;
                                    while (resize) {
                                        this_thread::yield();
                                    }
                                    cutWeights.resize(vectorSize);
                                    waitingForResize--;
                                }
                            }
                        }
                        result[Node] = bestCommunity;
#pragma omp atomic
                        communityVolumes[bestCommunity] += degree;
#pragma omp atomic
                        communityVolumes[currentCommunity] -= degree;
                        changed = true;
                        bool expected = true;
                        inQueue[Node].compare_exchange_strong(expected, false);
                        assert(expected);
                        graph.forNeighborsOf(Node, [&](node neighbor) {
                            if (result[neighbor] != bestCommunity && neighbor != Node) {      // Only add the node to the queue if it's not already in it, and it's not the Node we're currently moving
                                {
                                    expected = false;
                                    if (inQueue[neighbor].compare_exchange_strong(expected, true)) {
                                        newNodes.push_back(neighbor);
                                        if (newNodes.size() == WORKING_SIZE) {      //push new nodes to the queue once WORKING_SIZE has been reached
                                            qlock.lock();
                                            queue.emplace(move(newNodes));
                                            qlock.unlock();
                                            workAvailable.notify_all();         // Notify threads that new work is available
                                            newNodes.clear();
                                            newNodes.reserve(WORKING_SIZE);
                                        }
                                        assert(!expected);
                                    }
                                }
                            }
                        });
                    }

                }

                //queue check/wait
                // 3 cases : newnodes not empty -> continue, newnodes empty & queue not empty ->  pop queue, both empty -> increment waiting, if waiting < #threads wait, else done notify all
                PLIFF(totalNodesPerThread[omp_get_thread_num()] += currentNodes.size())
                if (!newNodes.empty()) {
                    swap(currentNodes, newNodes);
                    newNodes.clear();
                    continue;
                }

                unique_lock<mutex> uniqueLock(qlock);
                if (!queue.empty()) {
                    swap(currentNodes, queue.front());
                    queue.pop();
                } else {                                           // queue empty && newNodes empty
                    waitingForNodes++;
                    if (waitingForNodes < tcount) {                             // Not all nodes are done yet, wait for new work
                        waitingForResize++;
                        while (queue.empty() && waitingForNodes < tcount) {
                            workAvailable.wait(uniqueLock);
                        }
                        if (waitingForNodes < tcount) {                         // Notified and not all done means there's new work
                            swap(currentNodes, queue.front());
                            queue.pop();
                            waitingForNodes--;
                            waitingForResize--;
                            continue;
                        }
                    }
                    uniqueLock.unlock();                                        // Notified and all done, stop.
                    workAvailable.notify_all();
                    break;
                }
            } while (true);
#ifdef PLFINE
            #pragma omp critical(PL_TRACE)
            cout << "Thread " << omp_get_thread_num() << " worked " << totalNodesPerThread[omp_get_thread_num()] << "Nodes and moved " << moved[omp_get_thread_num()] << endl;
#endif
        }
        result.setUpperBound(upperBound);
        assert(queue.empty());
        assert(waitingForNodes == tcount);
#ifdef PLFINE
        long totalMoved = accumulate(moved.begin(), moved.end(), 0L);
        long totalWorked = accumulate(totalNodesPerThread.begin(), totalNodesPerThread.end(), 0L);
        cout << "Total worked: " << totalWorked << " Total moved: " << totalMoved << " moved to singleton community: " << singleton << endl;
#endif
    }

    Partition ParallelLeiden::parallelRefine(const Graph &graph) {
        Partition refined(graph.numberOfNodes());
        refined.allToSingletons();
#ifdef PLFINE
        cout << "Starting refinement with " << result.numberOfSubsets() << " partitions" << endl;
#endif
        vector<uint_fast8_t> singleton(refined.upperBound(), true);
        vector<double> cutCtoSminusC(refined.upperBound());
        vector<double> refinedVolumes(refined.upperBound());    // Community Volumes in the refined partition
        vector<mutex> locks(refined.upperBound());
        vector<node> nodes(graph.upperNodeIdBound(), none);
#pragma omp parallel
        {
            vector<index> neighComms;                             // Keeps track of relevant Neighbor communities. Needed to reset the clearlist fast
            vector<double> cutWeights(refined.upperBound());    // cut from Node to Communities
            auto &mt = Aux::Random::getURNG();
#pragma omp for
            for (node Node = 0; Node < graph.upperNodeIdBound(); Node++) {
                if (graph.hasNode(Node)) {
                    nodes[Node] = Node;
                    graph.forNeighborsOf(Node, [&](node neighbor, edgeweight ew) {
                        if (Node != neighbor) {
                            if (result[neighbor] == result[Node]) {             // Cut to communities in the refined partition that are in the same community in the original partition
                                cutCtoSminusC[Node] += ew;
                            }
                        } else {
                            refinedVolumes[Node] += ew;
                        }
                        refinedVolumes[Node] += ew;
                    });
                }
            }
            if (random) {
                int share = graph.upperNodeIdBound() / omp_get_num_threads();
                int start = omp_get_thread_num() * share;
                int end = (omp_get_thread_num() + 1) * share - 1;
                if (omp_get_thread_num() == omp_get_num_threads() - 1)
                    end = nodes.size() - 1;
                if (start != end && end > start)
                    shuffle(nodes.begin() + start, nodes.begin() + end, mt);
#pragma omp barrier
                }
#pragma omp for schedule(dynamic, WORKING_SIZE)
            for (node Node: nodes) {
                if (Node == none || !singleton[Node]) {              // only consider singletons
                    continue;
                }
                index S = result[Node];                               // Node's community ID in the previous partition (S)
                for (auto z: neighComms) {                             // Reset the clearlist : Set all cutweights to 0
                    if (z != none)
                        cutWeights[z] = 0;
                }

                neighComms.clear();

                vector<node> criticalNodes;                                     // Nodes whose community ID equals their Node ID. These are the only ones that can possibly
                double degree = 0;                                              // affect the cut which we need to update later since only those can be moved (possibly singletons)

                graph.forNeighborsOf(Node, [&](node neighbor, edgeweight ew) {      // Calculate degree and cut
                    degree += ew;
                    if (neighbor != Node) {
                        if (S == result[neighbor]) {
                            index z = refined[neighbor];
                            if (z == neighbor) {
                                criticalNodes.push_back(neighbor);      // Remember nodes that are still in their own Community. These might move and change the cut
                            }                                           // We don't need to remember the weight of that edge since it's already saved in cutWeights

                            if (cutWeights[z] == 0)
                                neighComms.push_back(z);        // Keep track of neighbor communities

                            cutWeights[z] += ew;
                        }
                    } else {
                        degree += ew;
                    }
                });
                if (cutCtoSminusC[Node] < this->gamma * degree * (communityVolumes[S] - degree) * inverseGraphVolume) {     // R-Set Condition
                    continue;
                }

                if (cutWeights[Node] != 0) { // Node has been moved -> not a singleton anymore. Stop.
                    continue;
                }

                double delta;
                index bestC = none;
                double bestDelta = numeric_limits<double>::lowest();
                int idx;
                auto bestCommunity = [&] {
                    for (int i = 0; i < neighComms.size(); i++) {         // Only consider (refined) communities the node is connected to
                        index C = neighComms[i];
                        if (C == none) {
                            continue;
                        }
                        delta = modularityDelta(cutWeights[C], degree, refinedVolumes[C]);

                        if (delta < 0) {              // modThreshold is 0, since cutw(v,C-) is 0 and volw(C-) is also 0
                            continue;
                        }

                        auto absC = refinedVolumes[C];
                        if (delta > bestDelta && cutCtoSminusC[C] >= this->gamma * absC * (communityVolumes[S] - absC) * inverseGraphVolume) { // T-Set Condition
                            bestDelta = delta;
                            bestC = C;
                            idx = i;
                        }
                    }
                };
                auto updateCut = [&] {
                    for (node &neighbor: criticalNodes) {
                        if (neighbor != none) {
                            index neighborCommunity = refined[neighbor];
                            if (neighborCommunity != neighbor) {
                                if (cutWeights[neighborCommunity] == 0) {
                                    neighComms.push_back(neighborCommunity);                      // remember to clear the vector, this community was not saved initially since the neighbor moved to it later
                                }
                                cutWeights[neighborCommunity] += cutWeights[neighbor];     // cutWeights[Neighbor] is the weight of the edge between Node and Neighbor, since Neighbor was a singleton
                                cutWeights[neighbor] = 0;                          // Clear cutWeights entry beforehand, so we can "erase" bestC from the pointers vector by replacing it with "none"
                                neighbor = none;
                            }
                        }
                    }
                };
                bestCommunity();
                if (bestC == none) {
                    continue;
                }
                lockLowerFirst(Node, bestC, locks);
                if (singleton[Node]) {                                  // If this node is no longer a singleton, stop.
                    while (bestC != none && refined[bestC] != bestC) {  // Target community still contains its "host" node? If not, then this community is now empty, choose a new one.
                        locks[bestC].unlock();
                        neighComms[idx] = none;                           // This makes sure it won't be considered in the next bestCommunity() call
                        bestC = none;
                        bestDelta = numeric_limits<double>::lowest();
                        updateCut();
                        bestCommunity();
                        if (bestC != none) {
                            if (!locks[bestC].try_lock()) {
                                if (Node < bestC) {
                                    locks[bestC].lock();
                                } else {
                                    locks[Node].unlock();                      // temporarily release lock on current Node to avoid deadlocks
                                    lockLowerFirst(Node, bestC, locks);
                                }
                                if (!singleton[Node]) {
                                    locks[Node].unlock();
                                    locks[bestC].unlock();
                                    continue;
                                }
                            }
                        }
                    }
                    if (bestC == none) {
                        locks[Node].unlock();           //bestC was already unlocked in the loop above
                        continue;
                    }
                    singleton[bestC] = false;
                    refined[Node] = bestC;
                    refinedVolumes[bestC] += degree;
                    updateCut();
                    cutCtoSminusC[bestC] += cutCtoSminusC[Node] - 2 * cutWeights[bestC];
                }
                locks[bestC].unlock();
                locks[Node].unlock();
            }
        }

#ifdef PLFINE
        cout << "Ending refinement with " << refined.numberOfSubsets() << " partitions" << endl;
#endif
        return refined;
    }
}