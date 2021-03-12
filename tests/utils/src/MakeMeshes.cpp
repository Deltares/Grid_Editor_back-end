#include <memory>
#include <stdexcept>

#include <boost/dll.hpp>

#include "../../../extern/netcdf/netCDF 4.6.1/include/netcdf.h"

#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>
#include <TestUtils/MakeMeshes.hpp>

Mesh2dPointers AllocateMesh2dData(meshkernelapi::Mesh2D mesh2d)
{
    std::unique_ptr<int> edge_nodes(new int[mesh2d.num_edges * 2]);
    std::unique_ptr<int> face_nodes(new int[mesh2d.num_face_nodes]);
    std::unique_ptr<int> nodes_per_face(new int[mesh2d.num_faces]);
    std::unique_ptr<double> node_x(new double[mesh2d.num_nodes]);
    std::unique_ptr<double> node_y(new double[mesh2d.num_nodes]);
    std::unique_ptr<double> edge_x(new double[mesh2d.num_edges]);
    std::unique_ptr<double> edge_y(new double[mesh2d.num_edges]);
    std::unique_ptr<double> face_x(new double[mesh2d.num_faces]);
    std::unique_ptr<double> face_y(new double[mesh2d.num_faces]);

    return {std::move(edge_nodes),
            std::move(face_nodes),
            std::move(nodes_per_face),
            std::move(node_x),
            std::move(node_y),
            std::move(edge_x),
            std::move(edge_y),
            std::move(face_x),
            std::move(face_y)};
}

meshkernelapi::Mesh2D ReadLegacyMeshFromFileForApiTesting(std::string filePath)
{

#if _WIN32
    boost::dll::shared_library lib("netcdf.dll");
#else
    boost::dll::shared_library lib("libnetcdf.so", boost::dll::load_mode::search_system_folders);
#endif

    auto nc_open = lib.get<int(const char*, int, int*)>("nc_open");
    auto nc_inq_dimid = lib.get<int(int, const char*, int*)>("nc_inq_dimid");
    auto nc_inq_dim = lib.get<int(int, int, char*, std::size_t*)>("nc_inq_dim");
    auto nc_inq_varid = lib.get<int(int, const char*, int*)>("nc_inq_varid");
    auto nc_get_var_double = lib.get<int(int, int, double*)>("nc_get_var_double");
    auto nc_get_var_int = lib.get<int(int, int, int*)>("nc_get_var_int");

    int ncidp = 0;
    int err = nc_open(filePath.c_str(), NC_NOWRITE, &ncidp);
    if (err != 0)
    {
        throw std::invalid_argument("ReadLegacyMeshFromFile: Could not load netcdf file.");
    }

    std::string mesh2dNodes{"nNetNode"};
    int dimid = 0;
    err = nc_inq_dimid(ncidp, mesh2dNodes.c_str(), &dimid);
    if (err != 0)
    {
        throw std::invalid_argument("ReadLegacyMeshFromFile: Could not find the ID of a dimension of 'nNetNode'.");
    }

    meshkernelapi::Mesh2D mesh2d{};

    std::size_t num_nodes;
    std::unique_ptr<char> read_name(new char[NC_MAX_NAME]);
    err = nc_inq_dim(ncidp, dimid, read_name.get(), &num_nodes);
    if (err != 0)
    {
        throw std::invalid_argument("ReadLegacyMeshFromFile: Could not gind the length of dimension of 'nNetNode'.");
    }
    mesh2d.num_nodes = static_cast<int>(num_nodes);

    std::string mesh2dEdges{"nNetLink"};
    err = nc_inq_dimid(ncidp, mesh2dEdges.c_str(), &dimid);
    if (err != 0)
    {
        throw std::invalid_argument("ReadLegacyMeshFromFile: Could not find the ID of a dimension of 'nNetLink'.");
    }

    std::size_t num_edges;
    err = nc_inq_dim(ncidp, dimid, read_name.get(), &num_edges);
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.node_x = new double[mesh2d.num_nodes];
    mesh2d.node_y = new double[mesh2d.num_nodes];

    std::string mesh2dNodeX{"NetNode_x"};
    int varid = 0;
    err = nc_inq_varid(ncidp, mesh2dNodeX.c_str(), &varid);
    err = nc_get_var_double(ncidp, varid, mesh2d.node_x);

    std::string mesh2dNodeY{"NetNode_y"};
    err = nc_inq_varid(ncidp, mesh2dNodeY.c_str(), &varid);
    err = nc_get_var_double(ncidp, varid, mesh2d.node_y);

    std::string mesh2dEdgeNodes{"NetLink"};
    err = nc_inq_varid(ncidp, mesh2dEdgeNodes.c_str(), &varid);

    mesh2d.edge_nodes = new int[mesh2d.num_edges * 2];
    err = nc_get_var_int(ncidp, varid, mesh2d.edge_nodes);

    // transform into 0 based indexing
    for (auto i = 0; i < mesh2d.num_edges * 2; i++)
    {
        mesh2d.edge_nodes[i] -= 1;
    }

    return mesh2d;
}

std::shared_ptr<meshkernel::Mesh2D> ReadLegacyMeshFromFile(std::string filePath, meshkernel::Projection projection)
{

    const auto mesh2d = ReadLegacyMeshFromFileForApiTesting(filePath);

    std::vector<meshkernel::Edge> edges(mesh2d.num_edges);
    std::vector<meshkernel::Point> nodes(mesh2d.num_nodes);

    for (auto i = 0; i < nodes.size(); i++)
    {
        nodes[i].x = mesh2d.node_x[i];
        nodes[i].y = mesh2d.node_y[i];
    }

    auto index = 0;
    for (auto i = 0; i < edges.size(); i++)
    {
        edges[i].first = mesh2d.edge_nodes[index];
        index++;
        edges[i].second = mesh2d.edge_nodes[index];
        index++;
    }

    auto mesh = std::make_shared<meshkernel::Mesh2D>(edges, nodes, projection);

    // clean up c memory
    DeleteRectangularMeshForApiTesting(mesh2d);

    return mesh;
}

std::shared_ptr<meshkernel::Mesh2D> MakeRectangularMeshForTesting(int n, int m, double delta, meshkernel::Projection projection, meshkernel::Point origin)
{
    std::vector<std::vector<size_t>> indicesValues(n, std::vector<size_t>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (auto i = 0; i < n; ++i)
    {
        for (auto j = 0; j < m; ++j)
        {
            indicesValues[i][j] = i * m + j;
            nodes[nodeIndex] = {origin.x + i * delta, origin.y + j * delta};
            nodeIndex++;
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;

    for (auto i = 0; i < n - 1; ++i)
    {
        for (auto j = 0; j < m; ++j)
        {
            edges[edgeIndex] = {indicesValues[i][j], indicesValues[i + 1][j]};
            edgeIndex++;
        }
    }

    for (auto i = 0; i < n; ++i)
    {
        for (auto j = 0; j < m - 1; ++j)
        {
            edges[edgeIndex] = {indicesValues[i][j + 1], indicesValues[i][j]};
            edgeIndex++;
        }
    }

    return std::make_shared<meshkernel::Mesh2D>(edges, nodes, projection);
}

meshkernelapi::Mesh2D MakeRectangularMeshForApiTesting(int n, int m, double delta)
{
    std::vector<std::vector<size_t>> indicesValues(n, std::vector<size_t>(m));
    meshkernelapi::Mesh2D mesh2d{};

    mesh2d.node_x = new double[n * m];
    mesh2d.node_y = new double[n * m];
    size_t nodeIndex = 0;
    for (auto i = 0; i < n; ++i)
    {
        for (auto j = 0; j < m; ++j)
        {

            mesh2d.node_x[nodeIndex] = i * delta;
            mesh2d.node_y[nodeIndex] = j * delta;
            indicesValues[i][j] = i * m + j;
            nodeIndex++;
        }
    }

    mesh2d.edge_nodes = new int[((n - 1) * m + (m - 1) * n) * 2];
    size_t edgeIndex = 0;
    for (auto i = 0; i < n - 1; ++i)
    {
        for (auto j = 0; j < m; ++j)
        {
            mesh2d.edge_nodes[edgeIndex] = indicesValues[i][j];
            edgeIndex++;
            mesh2d.edge_nodes[edgeIndex] = indicesValues[i + 1][j];
            edgeIndex++;
        }
    }

    for (auto i = 0; i < n; ++i)
    {
        for (auto j = 0; j < m - 1; ++j)
        {
            mesh2d.edge_nodes[edgeIndex] = indicesValues[i][j + 1];
            edgeIndex++;
            mesh2d.edge_nodes[edgeIndex] = indicesValues[i][j];
            edgeIndex++;
        }
    }

    mesh2d.num_nodes = nodeIndex;
    mesh2d.num_edges = edgeIndex / 2;

    return mesh2d;
}

void DeleteRectangularMeshForApiTesting(const meshkernelapi::Mesh2D& mesh2d)
{
    delete[] mesh2d.node_x;
    delete[] mesh2d.node_y;
    delete[] mesh2d.edge_nodes;
}

std::shared_ptr<meshkernel::Mesh2D> MakeSmallSizeTriangularMeshForTestingAsNcFile()
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({322.252624511719, 454.880187988281});
    nodes.push_back({227.002044677734, 360.379241943359});
    nodes.push_back({259.252227783203, 241.878051757813});
    nodes.push_back({428.003295898438, 210.377746582031});
    nodes.push_back({536.003967285156, 310.878753662109});
    nodes.push_back({503.753784179688, 432.379974365234});
    nodes.push_back({350.752807617188, 458.630249023438});
    nodes.push_back({343.15053976393, 406.232256102912});
    nodes.push_back({310.300984548069, 319.41005739802});
    nodes.push_back({423.569603308318, 326.17986967523});

    std::vector<meshkernel::Edge> edges;
    edges.push_back({2, 8});
    edges.push_back({1, 8});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    edges.push_back({3, 8});
    edges.push_back({1, 7});
    edges.push_back({0, 7});
    edges.push_back({0, 1});
    edges.push_back({7, 8});
    edges.push_back({6, 7});
    edges.push_back({0, 6});
    edges.push_back({8, 9});
    edges.push_back({7, 9});
    edges.push_back({3, 4});
    edges.push_back({4, 9});
    edges.push_back({3, 9});
    edges.push_back({5, 7});
    edges.push_back({5, 6});
    edges.push_back({5, 9});
    edges.push_back({4, 5});

    return std::make_shared<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);
}

std::shared_ptr<meshkernel::Mesh2D> MakeCurvilinearGridForTesting()
{
    std::vector<double> xCoordinates{777.2400642395231,
                                     776.8947176796199,
                                     776.5500495969297,
                                     776.2062753740686,
                                     775.8639152055595,
                                     789.0362557653892,
                                     788.6869233110746,
                                     788.3382299013459,
                                     787.9903970612199,
                                     787.6439255320029,
                                     800.8316506843166,
                                     800.4781249366497,
                                     800.125196925116,
                                     799.7730958077465,
                                     799.4223012497744,
                                     783.1381600024562,
                                     782.7908204953472,
                                     782.4441397491378,
                                     782.0983362176443,
                                     781.7539203687811,
                                     794.933953224853,
                                     794.5825241238622,
                                     794.231713413231,
                                     793.8817464344831,
                                     793.5331133908886,
                                     777.0673909595715,
                                     776.7223836382748,
                                     776.3781624854992,
                                     776.0350952898141,
                                     788.8615895382319,
                                     788.5125766062102,
                                     788.164313481283,
                                     787.8171612966114,
                                     800.6548878104832,
                                     800.3016609308829,
                                     799.9491463664313,
                                     799.5976985287605,
                                     782.9644902489017,
                                     782.6174801222426,
                                     782.2712379833911,
                                     781.9261282932127,
                                     794.7582386743576,
                                     794.4071187685466,
                                     794.0567299238571,
                                     793.7074299126859};

    std::vector<double> yCoordinates{1145.6125094106028,
                                     1147.65972166567,
                                     1149.7057403812394,
                                     1151.7507901691963,
                                     1153.7936194659985,
                                     1147.830596354202,
                                     1149.8579007347682,
                                     1151.884213062344,
                                     1153.9097092742677,
                                     1155.9332418178237,
                                     1150.0482786995756,
                                     1152.0565093840655,
                                     1154.063958611566,
                                     1156.0707463059343,
                                     1158.0758372662724,
                                     1146.7215528824024,
                                     1148.758811200219,
                                     1150.7949767217917,
                                     1152.830249721732,
                                     1154.8634306419112,
                                     1148.939437526889,
                                     1150.957205059417,
                                     1152.974085836955,
                                     1154.990227790101,
                                     1157.004539542048,
                                     1146.6361155381364,
                                     1148.6827310234548,
                                     1150.7282652752178,
                                     1152.7722048175974,
                                     1148.8442485444853,
                                     1150.8710568985562,
                                     1152.896961168306,
                                     1154.9214755460457,
                                     1151.0523940418207,
                                     1153.0602339978159,
                                     1155.06735245875,
                                     1157.0732917861033,
                                     1147.7401820413108,
                                     1149.7768939610055,
                                     1151.812613221762,
                                     1153.8468401818216,
                                     1149.948321293153,
                                     1151.965645448186,
                                     1153.982156813528,
                                     1155.9973836660745};

    std::vector<meshkernel::Edge> edges{{1, 16},
                                        {2, 17},
                                        {3, 18},
                                        {4, 19},
                                        {5, 20},
                                        {6, 21},
                                        {7, 22},
                                        {8, 23},
                                        {9, 24},
                                        {10, 25},
                                        {1, 26},
                                        {2, 27},
                                        {3, 28},
                                        {4, 29},
                                        {6, 30},
                                        {7, 31},
                                        {8, 32},
                                        {9, 33},
                                        {11, 34},
                                        {12, 35},
                                        {13, 36},
                                        {14, 37},
                                        {16, 38},
                                        {30, 38},
                                        {17, 38},
                                        {26, 38},
                                        {17, 39},
                                        {31, 39},
                                        {18, 39},
                                        {27, 39},
                                        {18, 40},
                                        {32, 40},
                                        {19, 40},
                                        {28, 40},
                                        {19, 41},
                                        {33, 41},
                                        {20, 41},
                                        {29, 41},
                                        {21, 42},
                                        {34, 42},
                                        {22, 42},
                                        {30, 42},
                                        {22, 43},
                                        {35, 43},
                                        {23, 43},
                                        {31, 43},
                                        {23, 44},
                                        {36, 44},
                                        {24, 44},
                                        {32, 44},
                                        {24, 45},
                                        {37, 45},
                                        {25, 45},
                                        {33, 45},
                                        {6, 16},
                                        {7, 17},
                                        {8, 18},
                                        {9, 19},
                                        {10, 20},
                                        {11, 21},
                                        {12, 22},
                                        {13, 23},
                                        {14, 24},
                                        {15, 25},
                                        {2, 26},
                                        {3, 27},
                                        {4, 28},
                                        {5, 29},
                                        {7, 30},
                                        {8, 31},
                                        {9, 32},
                                        {10, 33},
                                        {12, 34},
                                        {13, 35},
                                        {14, 36},
                                        {15, 37}};

    std::vector<meshkernel::Point> nodes(xCoordinates.size());

    for (auto i = 0; i < nodes.size(); i++)
    {
        nodes[i].x = xCoordinates[i];
        nodes[i].y = yCoordinates[i];
    }

    for (auto i = 0; i < edges.size(); i++)
    {
        edges[i].first -= 1;
        edges[i].second -= 1;
    }
    return std::make_shared<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);
}
