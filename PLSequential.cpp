#include "ParallelLeiden.hpp"

void ParallelLeiden::moveNodesFast(const Graph &graph) {
    queue<node> activeNodes;// keep track of active Nodes
    vector<node> newNodes;
    newNodes.reserve(graph.numberOfNodes() / 2);
    vector<bool> notInQueue(graph.numberOfNodes());   // Only insert nodes to the queue when they're not already in it.
    vector<double> cutWeights(result.upperBound() + 10000);  // cutWeight[Community] returns cut of Node to Community
    vector<index> pointers;

    graph.forNodesInRandomOrder([&](node Node) {
        activeNodes.push(Node);
    });

    do {
        node Node = activeNodes.front();
        activeNodes.pop();
        index currentCommunity = result[Node];                    // The node's current communityID
        double maxDelta = std::numeric_limits<double>::lowest();   // remember best delta so far ..
        index bestCommunity = none;                               // and the corresponding communityID
        double degree = 0;                                      // weighted degree of the Node
        notInQueue[Node] = true;
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

        for (auto community : pointers) {                                                         // For all neighbor communities, determine the modularity Delta for moving the node to the community
            if (community != currentCommunity) {                                                    // "Moving" a node to its current community is pointless
                double delta;
                delta = modularityDelta(cutWeights[community], degree, communityVolumes[community]);
                if (delta > maxDelta) {                                                             // Keep track of the best delta and the corresponding communityID
                    maxDelta = delta;
                    bestCommunity = community;
                }
            }
        }
        double modThreshold = modularityThreshold(cutWeights[currentCommunity], communityVolumes[currentCommunity], degree);
        if (0 > modThreshold || maxDelta > modThreshold) {
            if (0 > maxDelta) {
                result.toSingleton(Node);
                bestCommunity = result[Node];
                if (bestCommunity > communityVolumes.capacity()) {
                    communityVolumes.resize(bestCommunity + 10000);
                    cutWeights.resize(bestCommunity + 10000);
                }
            } else {
                result[Node] = bestCommunity;
            }
            communityVolumes[bestCommunity] += degree;
            communityVolumes[currentCommunity] -= degree;
            changed = true;
            graph.forNeighborsOf(Node, [&](node neighbor) {
                if (notInQueue[neighbor] && result[neighbor] != bestCommunity && neighbor != Node) {      // Only add the node to the queue if it's not already in it, and it's not the Node we're currently moving
                    {
                        notInQueue[neighbor] = false;
                        activeNodes.push(neighbor);
                    }
                }
            });


        }
        for (auto z : pointers) {                             // Reset the clearlist : Set all cutweights to 0 and clear the pointer vector
            cutWeights[z] = 0;
        }
        pointers.clear();
    } while (!activeNodes.empty());
}

Partition ParallelLeiden::refineAndMerge(const Graph &graph) {
    Partition refined(graph.numberOfNodes());
    refined.allToSingletons();
#ifdef PLDEBUG
    cout << "Starting refinement with " << result.numberOfSubsets() << " partitions" << endl;
#endif
    vector<bool> singleton(refined.upperBound(), true);
    vector<double> cutWeights(refined.upperBound());
    vector<double> cutCtoSminusC(refined.upperBound());
    vector<double> refinedVolumes(refined.upperBound());
    vector<index> pointers;
    auto &mt = Aux::Random::getURNG();

    graph.forNodes([&](node Node) {
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
    });

    graph.forNodesInRandomOrder([&](node Node) {
        if (!singleton[refined[Node]]) {              // only consider singletons
            return;
        }
        index S = result[Node];                                                 // Node's community ID in the previous partition (S)
        index refinedCommunity = refined[Node];                                   // Node's community ID in the refined partition (refined)

        for (auto z : pointers) {                             // Reset the clearlist : Set all cutweights to 0
            cutWeights[z] = 0;
        }
        pointers.clear();

        double degree = 0;
        graph.forNeighborsOf(Node, [&](node neighbor, edgeweight ew) {      // Calculate degree and cut to S-v
            index neighborCommunity = result[neighbor];
            degree += ew;
            if (neighbor != Node) {
                if (S == neighborCommunity) {
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
            return;
        }

        double delta;
        double modularityThreshold = cutWeights[refinedCommunity] - (refinedVolumes[refinedCommunity] - degree) * gamma * degree * inverseGraphVolume;
        index bestC = none;
        double bestDelta = 0;
        for (const auto C : pointers) {         // Only consider (refined) communities the node is connected to
            delta = modularityDelta(cutWeights[C], degree, refinedVolumes[C]);

            if (delta < modularityThreshold) {              // actual modularity delta >= 0
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
            return;
        }
        singleton[bestC] = false;
        refinedVolumes[bestC] += degree;
        cutCtoSminusC[bestC] += cutCtoSminusC[refinedCommunity] - 2 * cutWeights[bestC];
        refined.moveToSubset(bestC, Node);


    });
#ifdef PLDEBUG
    cout << "Ending refinement with " << refined.numberOfSubsets() << " partitions" << endl;
#endif
    return refined;
}