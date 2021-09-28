#include "ParallelLeiden.hpp"
#define WORKING_SIZE 1000

/**
 *      As the filename already suggest these functions are EXPERIMENTAL,
 *      they have NOT been tested thoroughly and should NOT be used unless you know what you're doing
 *      In particular the resizing of the communityVolume and cutWeights vectors has been blindly
 *      transferred from the parallelMove function. These functions may deadlock or simply return wrong data.
 */


void NetworKit::ParallelLeiden::experimentalMove(const Graph &graph) {
#ifdef PLDEBUG
    cout << "Local Moving : " << graph.numberOfNodes() << " Nodes " << endl;
#endif
    vector<atomic_bool> inQueue(graph.upperNodeIdBound());   // Only insert nodes to the queue when they're not already in it
#pragma omp parallel
    {
        vector<node> newNodes;
        vector<node> activeNodes;
        activeNodes.reserve(graph.numberOfNodes() / omp_get_num_threads() + 1);
        vector<double> cutWeights(communityVolumes.capacity());  // cutWeight[Community] returns cut of Node to Community
        vector<index> pointers;
        atomic_bool resize = false;
        atomic_int waitingForResize = 0;
        uint64_t vectorSize = communityVolumes.capacity();
#pragma omp for
        for (node Node = 0; Node < graph.upperNodeIdBound(); Node++) {
            if (graph.hasNode(Node)) {
                activeNodes.push_back(Node);
                atomic_init(&inQueue[Node], true);
            }
        }
        if (random)
            shuffle(activeNodes.begin(), activeNodes.end(), Aux::Random::getURNG());
#pragma omp barrier
        do {
            for (node Node: activeNodes) {
                if (resize) {                                                   // Likely wont happen, like, EVER.
                    waitingForResize++;
                    while (resize) {
                        this_thread::yield();
                    }
                    waitingForResize--;
                }
                cutWeights.resize(vectorSize);
                for (auto z: pointers) {                             // Reset the clearlist : Set all cutweights to 0 and clear the pointer vector
                    cutWeights[z] = 0;
                }
                pointers.clear();
                index currentCommunity = result[Node];                    // The node's current communityID
                double maxDelta = -std::numeric_limits<double>::max();   // remember best delta so far ..
                index bestCommunity = none;                               // and the corresponding communityID
                double degree = 0;                                      // weighted degree of the Node
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
                    degree += ew;                                       // keep track of the nodes degree, since we're already iterating over its neighbors.

                });

                double delta;
                for (auto community: pointers) {                                                         // For all neighbor communities, determine the modularity Delta for moving the node to the community
                    if (community != currentCommunity) {                                                    // "Moving" a node to its current community is pointless
                        delta = modularityDelta(cutWeights[community], degree, communityVolumes[community]);
                        if (delta > maxDelta) {                                                             // Keep track of the best delta and the corresponding communityID
                            maxDelta = delta;
                            bestCommunity = community;
                        }
                    }
                }
                double modThreshold = modularityThreshold(cutWeights[currentCommunity], communityVolumes[currentCommunity], degree);
                if (0 > modThreshold || maxDelta > modThreshold) {
                    if (0 > maxDelta) {                         // move node to empty community
                        bestCommunity = result.upperBound();
                        if (bestCommunity >= communityVolumes.capacity()) {                 // Chances are this will never happen, EVER. Even Graphs with 100million nodes aren't close to hitting this.
                            bool expected = false;
                            if (resize.compare_exchange_strong(expected, true)) {
                                vectorSize += 10000;
                                cutWeights.resize(vectorSize);
                                while (waitingForResize < omp_get_num_threads() - 1) {
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
#pragma omp critical (PL_singleton)                                        // This happens extremely few times and only towards the end of an iteration of the algorithm
                        result.toSingleton(Node);
                        bestCommunity = result[Node];
                    } else {
                        result[Node] = bestCommunity;
                    }
#pragma omp atomic
                    communityVolumes[bestCommunity] += degree;
#pragma omp atomic
                    communityVolumes[currentCommunity] -= degree;
                    bool expected = true;
                    inQueue[Node].compare_exchange_strong(expected, false);
                    assert(expected);
                    changed = true;
                    graph.forNeighborsOf(Node, [&](node neighbor) {
                        if (result[neighbor] != bestCommunity && neighbor != Node) {      // Only add the node to the queue if it's not already in it, and it's not the Node we're currently moving
                            {
                                expected = false;
                                if (inQueue[neighbor].compare_exchange_strong(expected, true)) {
                                    newNodes.push_back(neighbor);
                                }
                            }
                        }
                    });
                }
            }
            activeNodes.clear();
            swap(activeNodes, newNodes);
        } while (!activeNodes.empty());
        waitingForResize++;
    }
}

Partition ParallelLeiden::experimentalRefine(const Graph &graph) {
    PLPRINT("Starting refinement with " << result.numberOfSubsets() << " partitions")
    Partition refined(graph.numberOfNodes());
    refined.allToSingletons();
    vector<uint8_t> singleton(graph.upperNodeIdBound(), true);         // Avoid specialized bool vector
    vector<double> cutCtoSminusC(refined.upperBound());
    vector<double> refinedVolumes(graph.upperNodeIdBound());

    vector<pair<node, index>> nodes;
    nodes.reserve(graph.numberOfNodes());
    graph.forNodes([&](node Node) {
        nodes.emplace_back(Node, result[Node]);
    });
    auto cmp = [&](const pair<node, index> &p1, const pair<node, index> &p2) {
        return p1.second < p2.second;
    };
    auto cmpsize = [&](const pair<node, index> &p1, const pair<node, index> &p2) {
        return p1.second - p1.first > p2.second - p2.first;
    };
    Aux::Parallel::sort(nodes.begin(), nodes.end(), cmp);
    vector<pair<index, index>> intervals;
    index prev = 0;
    index idx = nodes[0].second;
    for (int i = 1; i < graph.numberOfNodes(); i++) {
        if (nodes[i].second != idx) {
            intervals.emplace_back(prev, i);
            prev = i;
            idx = nodes[i].second;
        }
    }
    intervals.emplace_back(prev, nodes.size());
    if (random)
        shuffle(intervals.begin(), intervals.end(), Aux::Random::getURNG());
    else
        Aux::Parallel::sort(intervals.begin(), intervals.end(), cmpsize);
    atomic_int crt;
#pragma omp parallel
    {
#pragma omp single
        crt = omp_get_num_threads();
#pragma omp for
        for (node Node = 0; Node < graph.upperNodeIdBound(); Node++) {
            if (!graph.hasNode(Node))
                continue;
            graph.forNeighborsOf(Node, [&](node neighbor, edgeweight ew) {
                if (Node != neighbor) {
                    if (result[neighbor] == result[Node]) {                             // Cut to communities in the refined partition that are in the same community in the original partition
                        cutCtoSminusC[Node] += ew;
                    }
                } else {
                    refinedVolumes[Node] += ew;
                }
                refinedVolumes[Node] += ew;
            });
        }

        vector<double> cutWeights(refined.upperBound());
        vector<index> pointers;
        int z = omp_get_thread_num();
        while (z < intervals.size()) {
            for (int i = intervals[z].first; i < intervals[z].second; i++) {
                node Node = nodes[i].first;
                if (!singleton[Node]) {              // only consider singletons
                    continue;
                }
                index S = result[Node];                                                 // Node's community ID in the previous partition (S)

                for (auto z: pointers) {                             // Reset the clearlist : Set all cutweights to 0
                    cutWeights[z] = 0;
                }
                pointers.clear();
                double degree = 0;
                graph.forNeighborsOf(Node, [&](node neighbor, edgeweight ew) {      // Calculate degree and cut to S-v
                    index neighborCluster = result[neighbor];
                    degree += ew;
                    if (neighbor != Node) {
                        if (S == neighborCluster) {

                            index z = refined[neighbor];
                            if (cutWeights[z] == 0) {
                                pointers.push_back(z);        // Keep track of neighbor communities
                            }
                            cutWeights[z] += ew;
                        }
                    } else {
                        degree += ew;
                    }
                });
                if (cutCtoSminusC[Node] < this->gamma * degree * (communityVolumes[S] - degree) * inverseGraphVolume) {     // R-Set Condition
                    continue;
                }

                double delta;
                index bestC = none;
                double bestDelta = 0;
                for (const auto C: pointers) {         // Only consider (refined) communities the node is connected to
                    delta = modularityDelta(cutWeights[C], degree, refinedVolumes[C]);

                    if (delta < 0) {              // actual modularity delta >= 0
                        continue;
                    }

                    auto absC = refinedVolumes[C];
                    if (cutCtoSminusC[C] >= this->gamma * absC * (communityVolumes[S] - absC) * inverseGraphVolume) {
                        if (delta > bestDelta) {
                            bestDelta = delta;
                            bestC = C;
                        }
                    }
                }

                if (bestC == none) {
                    continue;
                }
                singleton[bestC] = false;
                refinedVolumes[bestC] += degree;
                cutCtoSminusC[bestC] += cutCtoSminusC[Node] - 2 * cutWeights[bestC];
                refined[Node] = bestC;
            }
            z = crt++;
        }
    }
#ifdef PLDEBUG
    cout << "Ending refinement with " << refined.numberOfSubsets() << " partitions" << endl;
#endif
    return refined;
}
