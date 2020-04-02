//
//  main.cpp
//  bitSetClusteringProblem
//
//  Created by Utkarsh Khandelwal on 02/04/20.
//  Copyright © 2020 Utk_Sort. All rights reserved.
//
//
//  bitset_Problem.cpp
//  ClusteringAlgorithm
//
//  Created by Utkarsh Khandelwal on 02/04/20.
//  Copyright © 2020 Utk_Sort. All rights reserved.
//

#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <array>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>
#include <queue>
#include <bitset>
#include <chrono>

using integer = long long int;
using uinteger = unsigned long long int;

using namespace std;


static constexpr uinteger bitsCount = 24;
struct Edge;
struct Node
{
   integer m_id;
   shared_ptr<Node> m_leader;
   uinteger m_clusterSize;
   std::vector<shared_ptr<Edge>> m_allEdges;
   std::string m_bitStr;
   std::bitset<bitsCount> m_bits;
};

struct Edge
{
    shared_ptr<Node> e_nodeA;
    shared_ptr<Node> e_nodeB;
    integer e_cost;
};

uinteger InterfaceBits(vector<shared_ptr<Node>> &nodes, vector<shared_ptr<Edge>> &edges, string &inputFilePath)
{
    //Hash Table to check for nodes already added;


    
    uinteger nodesCount  = 0;
    integer bitsUsed = 0;

    ifstream input(inputFilePath);
    string lineData;
    uinteger linesRead = 0;
    while(getline(input, lineData))
    {
        if(linesRead == 0)
        {
            std::vector<integer> lineDataValues;
            stringstream strm(lineData);
            integer value;
            while(strm >> value)
            {
                lineDataValues.push_back(value);
            }
            nodesCount = lineDataValues[0];
            bitsUsed = lineDataValues[1];
            linesRead++;
        }
        else
        {
            linesRead++;
            std::string bitStr;
            
            stringstream sstream(lineData);
            std::string temp;
            while(sstream >> temp)
            {
                bitStr.append(temp);
            }
            
            
            std::bitset<bitsCount> bitPattern(bitStr);
            
            shared_ptr<Node> newNode = make_shared<Node>();
            newNode->m_id = linesRead;
            newNode->m_leader = newNode;
            newNode->m_clusterSize = 1;
            newNode->m_bitStr = bitStr;
            newNode->m_bits = bitPattern;
            nodes.push_back(newNode);
        }
    }
    
    //Now creating edges
    for(uinteger iNode = 0; iNode < nodesCount; ++iNode)
    {
        for(uinteger jNode = iNode + 1; jNode < nodesCount; ++jNode)
        {
            shared_ptr<Node> node1 = nodes[iNode];
            shared_ptr<Node> node2 = nodes[jNode];
            
            auto xor_result = node1->m_bits ^ node2->m_bits;
            if(xor_result.count() <=2)
            {
                //Create and store that edge
                std::shared_ptr<Edge> currEdge = make_shared<Edge>();
                currEdge->e_cost = xor_result.count();
                currEdge->e_nodeA = node1;
                currEdge->e_nodeB = node2;
                edges.push_back(currEdge);
            }
        }
    }
    return nodesCount;
}

bool FindConnetedNodes(vector<shared_ptr<Node>> &connectedNodes, shared_ptr<Node> &seedNode)  //BFS
{
    std::unordered_set<shared_ptr<Node>> visitedNodes;
    std::queue<shared_ptr<Node>> queue;
    queue.push(seedNode);
    visitedNodes.insert(seedNode);
    connectedNodes.push_back(seedNode);

    while(queue.size() != 0)
    {
        auto node = queue.front();
        auto &edges = node->m_allEdges;
        for(auto edge : edges)
        {
            shared_ptr<Node> otherNode;
            if(edge->e_nodeA == node)
                otherNode = edge->e_nodeB;
            else
                otherNode = edge->e_nodeA;
            if(visitedNodes.find(otherNode) != visitedNodes.end())
            {
                continue;
            }
            else
            {
                visitedNodes.insert(otherNode);
                connectedNodes.push_back(otherNode);
                queue.push(otherNode);
            }
        }
        queue.pop();
    }
    return true;
}


bool JoinClusters(shared_ptr<Edge> &currEdge, shared_ptr<Node> &node1, shared_ptr<Node> &node2)
{
    uinteger newClusterSize = node1->m_clusterSize + node2->m_clusterSize;
    std::unordered_set<shared_ptr<Node>> markedNodes;

    shared_ptr<Node> clusterToUpdate, leaderNode;

    if(node1->m_leader->m_clusterSize >= node2->m_leader->m_clusterSize)
    {
        clusterToUpdate = node2;
        leaderNode = node1->m_leader;
    }
    else
    {
        clusterToUpdate = node1;
        leaderNode = node2->m_leader;
    }

    std::vector<shared_ptr<Node>> connectedNodes;
    FindConnetedNodes(connectedNodes, clusterToUpdate);
    for(auto node : connectedNodes)
    {
        node->m_leader = leaderNode;
    }
    leaderNode->m_clusterSize = newClusterSize;
    node2->m_allEdges.push_back(currEdge);
    node1->m_allEdges.push_back(currEdge);
    return true;
}

bool CheckIfBelongsToSameCluster(shared_ptr<Node> &node1, shared_ptr<Node> &node2)
{
    if(node1->m_leader->m_id == node2->m_leader->m_id)
        return true;
    else
        return false;
}

integer ReduceClusters(uinteger currClusterCount, std::vector<shared_ptr<Edge>> &edges, vector<shared_ptr<Node>> &allNodes)
{
    auto distances = edges;
    std::sort(distances.begin(), distances.end(), [](const shared_ptr<Edge> &edge1, const shared_ptr<Edge> &edge2)->bool
                                                    {
                                                        if(edge1->e_cost < edge2->e_cost)
                                                            return true;
                                                        else
                                                            return false;
                                                    });

//    //For debugging store distances
//    vector<integer> edgeCosts;
//    for(auto dis : distances)
//    {
//        edgeCosts.push_back(dis->e_cost);
//    }

    uinteger iEdge = 0, lastEdgeUsed = 0;
    for(; iEdge < distances.size(); ++iEdge)
    {
        auto currEdge = distances[iEdge];
        auto node1 = currEdge->e_nodeA;
        auto node2 = currEdge->e_nodeB;
        bool belongsToSameCluster = CheckIfBelongsToSameCluster(node1, node2);

        if(!belongsToSameCluster)
        {
             JoinClusters(currEdge, node1, node2);
            currClusterCount--;
        }

//        //For debugging
//        cout << "iEdge:  " << iEdge << std::endl;
//        for(auto currNode : allNodes)
//        {
//            cout << "Node Value:  " << currNode->m_id << "   Leader Value:   " << currNode->m_leader->m_id << std::endl;
//        }
    }
    return currClusterCount;
}



int main(int argc, const char * argv[]) {


    vector<shared_ptr<Edge>> graphEdges;
    vector<shared_ptr<Node>> graphNodes;
    //string filePath("/Users/utkarsh/Algorithms_Learning/Sorting_Alogs_CPP/Sorting_Algorithms_Learning/stanford-algs-master/testCases/course3/assignment2Clustering/question2/input_random_72_16384_24.txt");

    string filePath("/Users/utkarsh/Downloads/clustering_big.txt");
    auto timePoint1 = std::chrono::high_resolution_clock::now();
    uinteger clustersCount = InterfaceBits(graphNodes, graphEdges, filePath);
    auto timePoint2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> creatingGraph = timePoint2 - timePoint1;
    std::cout << "Time taken to load :    " << creatingGraph.count() << "milliSec" << std::endl;
    //Just for debugging
    //edges counts
    
    uinteger graphEdgesCount = graphEdges.size();
    
    auto timePoint3 = std::chrono::high_resolution_clock::now();
    integer remainingClusterCounts = ReduceClusters(clustersCount, graphEdges, graphNodes);
    auto timePoint4 = std::chrono::high_resolution_clock::now();



    std::chrono::duration<double, std::milli> clustering = timePoint4 - timePoint3;
    std::cout << "Number of Clusters :    " << remainingClusterCounts << std::endl;

    std::cout << "Time taken to load :    " << clustering.count() << "milliSec" << std::endl;
    
    return 0;
}


