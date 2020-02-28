#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../Polygons.hpp"
#include "../Constants.cpp"
#include "../Orthogonalization.cpp"
#include "../GeometryListNative.hpp"
#include "../Splines.hpp"
#include <gtest/gtest.h>
#include <chrono>

#if defined(_WIN32)
#include <Windows.h>
#endif

TEST(OrthogonalizationTests, TestOrthogonalizationOneQuadOneTriangle)
{
    // Preparation
    std::vector<GridGeom::Point> nodes;
    nodes.push_back(GridGeom::Point{ 0.0,0.0 });
    nodes.push_back(GridGeom::Point{ 0.0,10.0 });
    nodes.push_back(GridGeom::Point{ 10.0,0.0 });
    nodes.push_back(GridGeom::Point{ 10.0,10.0 });
    nodes.push_back(GridGeom::Point{ 20.0,0.0 });

    std::vector<GridGeom::Edge> edges;
    edges.push_back({ 1, 0 });
    edges.push_back({ 0, 2 });
    edges.push_back({ 2, 4 });
    edges.push_back({ 4, 3 });
    edges.push_back({ 3, 2 });
    edges.push_back({ 3, 1 });

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 0;
    int projectToLandBoundaryOption = 0;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    GridGeomApi::GeometryListNative geometryListNativePolygon;
    GridGeomApi::GeometryListNative geometryListNativeLandBoundaries;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;

    // Execute
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    GridGeom::Orthogonalization orthogonalization;
    std::vector<GridGeom::Point> polygon;
    std::vector<GridGeom::Point> landBoundary;

    orthogonalization.Set(mesh, 
        isTriangulationRequired, 
        isAccountingForLandBoundariesRequired, 
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);

    orthogonalization.Iterate(mesh);

    // Assert
    constexpr double tolerance = 1e-8;
    ASSERT_NEAR(0.0, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(20.0, mesh.m_nodes[4].x, tolerance);

    ASSERT_NEAR(0.0, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[4].y, tolerance);
}

TEST(OrthogonalizationTests, TestOrthogonalizationSmallTriangularGrid)
{
    //One gets the edges
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 322.252624511719,454.880187988281 });
    nodes.push_back(GridGeom::Point{ 227.002044677734,360.379241943359 });
    nodes.push_back(GridGeom::Point{ 259.252227783203,241.878051757813 });
    nodes.push_back(GridGeom::Point{ 428.003295898438,210.377746582031 });
    nodes.push_back(GridGeom::Point{ 536.003967285156,310.878753662109 });
    nodes.push_back(GridGeom::Point{ 503.753784179688,432.379974365234 });
    nodes.push_back(GridGeom::Point{ 350.752807617188,458.630249023438 });
    nodes.push_back(GridGeom::Point{ 343.15053976393,406.232256102912 });
    nodes.push_back(GridGeom::Point{ 310.300984548069,319.41005739802 });
    nodes.push_back(GridGeom::Point{ 423.569603308318,326.17986967523 });

    std::vector<GridGeom::Edge> edges;
    // Local edges
    edges.push_back({ 3, 9 });
    edges.push_back({ 9, 2 });
    edges.push_back({ 2, 3 });
    edges.push_back({ 3, 4 });
    edges.push_back({ 4, 9 });
    edges.push_back({ 2, 8 });
    edges.push_back({ 8, 1 });
    edges.push_back({ 1, 2 });
    edges.push_back({ 9, 8 });
    edges.push_back({ 8, 7 });
    edges.push_back({ 7, 1 });
    edges.push_back({ 9, 10 });
    edges.push_back({ 10, 8 });
    edges.push_back({ 4, 5 });
    edges.push_back({ 5, 10 });
    edges.push_back({ 10, 4 });
    edges.push_back({ 8, 6 });
    edges.push_back({ 6, 7 });
    edges.push_back({ 10, 6 });
    edges.push_back({ 5, 6 });

    for (int i = 0; i < edges.size(); i++)
    {
        edges[i].first -= 1;
        edges[i].second -= 1;
    }

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 0;
    int projectToLandBoundaryOption = 0;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    GridGeomApi::GeometryListNative geometryListNativePolygon;
    GridGeomApi::GeometryListNative geometryListNativeLandBoundaries;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;

    // now build node-edge mapping
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);
    GridGeom::Orthogonalization orthogonalization;
    std::vector<GridGeom::Point> polygon;
    std::vector<GridGeom::Point> landBoundary;

    orthogonalization.Set(mesh,
        isTriangulationRequired,
        isAccountingForLandBoundariesRequired,
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);

    orthogonalization.Iterate(mesh);

    constexpr double tolerance = 1e-2;

    ASSERT_NEAR(325.590101919525, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(229.213730481198, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(263.439319753147, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(429.191105834504, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(535.865215426468, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(503.753784179688, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(354.048340705929, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(346.790050854504, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(315.030130405285, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(424.314957449766, mesh.m_nodes[9].x, tolerance);

    ASSERT_NEAR(455.319334078551, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(362.573521507281, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(241.096458631763, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(211.483073921775, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(311.401495506714, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(432.379974365234, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(458.064836627594, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(405.311585650679, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(319.612138503550, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(327.102805172725, mesh.m_nodes[9].y, tolerance);
}

TEST(OrthogonalizationTests, TestOrthogonalizationMediumTriangularGrid)
{
    std::vector<double> xCoordinates{ 62.7625648300453, 161.949072158728, 261.069585881221,
    360.172666249972, 459.142369029805, 558.084265626781, 656.824528757364,
    755.528481353452, 853.965329311203, 952.359934686054, 1050.43160060353,
    1120.05712659647, 1189.11500645678, 1258.14163387914, 1326.90596004565,
    1395.64219936137, 1464.13134631492, 1532.59797942102, 1600.85544140631,
    1669.09790464577, 1668.80002792144, 1576.70387218423, 1484.72477563802,
    1392.74180604902, 1300.78623773346, 1208.8152066093, 1116.79546197923,
    1024.75051485453, 932.604068150081, 840.425372250284, 748.1197979515,
    655.77810127546, 563.309505067612, 493.496955812432, 440.669090313994,
    387.821668971785, 334.92821561367, 282.018066250363, 229.078085651579,
    176.12540276634, 123.156212750692, 70.179099233301, 17.1950050344739,
    -36.4298784330279, 1602.78943256412, 1493.57352970632, 571.84959799321,
    406.808939831148, 137.797902110879, 1003.3765605098, 1023.76408459998,
    1023.35696842049, 1517.00160392216, 1449.83685833077, 1380.82852308108,
    1549.75220853121, 1466.66788680045, 505.44768153021, 517.720727684875,
    488.900806728985, 595.534139457637, 640.283702687422, 613.159414848467,
    555.316577425639, 300.985129199092, 296.284993710995, 325.161864449324,
    123.721554470859, 73.7435999462389, 208.083637673065, 218.942694781924,
    305.06983243846, 229.493547515462, 195.79033801243, 158.078345459855,
    264.2826875886, 210.446742542291, 305.63163904094, 438.190917250183,
    365.384481615292, 376.108055692517, 408.542989862293, 364.946483616079,
    1121.06618536549, 819.020359815983, 936.930762199942, 1039.75557546871,
    1094.37872969698, 848.036134203893, 1162.25791777307, 1008.49493002575,
    925.813108723684, 1330.6341283351, 1300.00925949891, 1440.50752293259,
    1392.61948009526, 1528.61628200069, 1454.92113805161, 1398.74369926468,
    1341.61181164575, 1386.56283193265, 1305.37214075439, 1241.3157231037,
    1308.75366968582, 1256.48584317216, 1168.29404203633, 1120.66341388579,
    1241.75084879338, 1214.55123478042, 1195.48628034708, 1190.54189328242,
    1144.32019801741, 1036.47382051833, 1066.38529141324, 1099.52012660602,
    1104.90451393637, 923.139518803558, 992.387479354801, 914.797632142623,
    1005.55684894905, 969.219626178607, 1020.3245941422, 1005.58030633185,
    936.786291497255, 967.496370997864, 900.506065618184, 848.992313167665,
    880.73960815057, 766.04287959682, 836.64987962007, 738.986567970715,
    884.330967563659, 808.148320217935, 745.443247488341, 660.796052447892,
    658.224804896114, 825.555372121912, 796.559766313242, 727.545484738486,
    860.977282078688, 914.732737725336, 440.384401490909, 556.117614932609,
    717.547231468049, 530.107250307386, 239.724241284403, 232.369957858814,
    277.281659520074, 344.550151745036, 320.693511828882, 378.604778671573,
    466.318457755492, 466.122785330259, 127.874606183712, 313.695801434773,
    1146.2620934797, 1196.72614909897, 1263.04312192961, 1331.93205293064,
    855.452236434789, 962.141515608211, 1058.82204784257, 1083.09776014187,
    984.275892728239, 894.710252045683, 1135.78919862608, 1078.07874733901,
    1035.96418590343, 733.783275720045, 1232.5723916232, 966.544577257275,
    910.471604720219, 1120.90869472635, 1170.05842291937, 1105.53580143792,
    1169.53766862171, 607.573206506069, 521.086178873613, 589.484624647177,
    525.288406830452, 650.444502107815, 479.132558194633, 960.382625206505,
    389.381254965268, 444.210045282223, 385.555611173953, 468.715850599873,
    659.895534601595, 700.896923431698, 672.750112485226, 749.608966218665,
    596.493635228729, 786.522274798097, 684.442960165764, 649.702994733637,
    644.7183213309, 763.369879515855, 560.913867115511, 706.979155959597,
    751.94594530305, 707.744280870667, 798.329487341725, 715.607133439597,
    626.150308205423, 515.323772235382, 608.758850892345, 683.017850450226,
    661.500134676386, 587.883140228208, 404.65011295109, 874.56975518804,
    1134.81205618378, 974.008894493711, 800.503884711923, 716.852262243643,
    1228.89709618001, 1289.21663648852, 1285.31338191788, 1168.37560582403,
    1071.73988023703, 1052.50389873111, 494.910294025769, 468.424842667466,
    501.1653800145, 706.697699611797, 587.48668380637, 567.591166719208,
    352.955968846849, 519.183613864093, 545.549863431857, 475.473855647218,
    575.93266335685, 779.63478167254, 783.010981800924, 856.079776531503,
    881.687564699277, 924.808161528524, 686.591364483249, 582.423948290607,
    1016.82100408926, 947.440193183191, 1102.67075930161, 442.397615031199,
    496.553903344142, 816.391475157116, 781.091077132163, 824.459666141569,
    405.386251631545 };

    std::vector<double> yCoordinates{ 1398.17508021161, 1424.58568867191, 1450.97872487021,
    1477.36711905124, 1503.71999848965, 1530.06547391488, 1556.35725995557,
    1582.63937751074, 1608.85037252908, 1635.05011952238, 1661.16387671934,
    1588.18582878269, 1514.38868265746, 1440.62493374534, 1367.14148715564,
    1293.6880549457, 1220.4986724169, 1147.33334880353, 1074.39155119179,
    1001.46578165158, 943.525148452222, 905.485431952072, 867.494065987331,
    829.501100287529, 791.519452505016, 753.531417910251, 715.523262519572,
    677.504697402845, 639.444208546661, 601.370399370658, 563.244183899422,
    525.10304831584, 486.90949770825, 519.084483046858, 604.230572144341,
    689.408180660606, 774.659981955451, 859.938693282193, 945.265485541408,
    1030.61275089762, 1115.98662186402, 1201.37326365087, 1286.77115665368,
    1371.76289117706, 972.83559570115, 976.409566518165, 597.265927654306,
    996.450070651327, 1302.46215356866, 859.684279104672, 1428.86668547679,
    1565.63048854101, 919.688965785097, 925.145712647591, 1028.31069245344,
    1022.62196875989, 1066.3710986491, 585.460417885342, 709.075452600758,
    649.520043085345, 539.924406214115, 591.970054875036, 686.271814751933,
    1341.49604828113, 1167.19944785855, 1046.06376207114, 945.793126071363,
    1368.69573155234, 1323.46154126497, 1352.01227233799, 1409.86578868729,
    1363.75392008515, 1247.01836238194, 1294.82561147764, 1232.45567928892,
    1307.32773866493, 1158.03270103058, 1243.28716538548, 1197.00632648838,
    1201.45245564987, 1103.62560745703, 1309.84889850827, 1264.41694036462,
    835.580016568779, 783.814702783565, 759.544997601521, 771.608769204143,
    1494.94644510885, 1433.16977103559, 1281.9830299247, 1621.44569940126,
    1536.88512232497, 1164.21474221973, 923.706016925083, 1009.54284406983,
    961.535054106462, 1083.25017619028, 1143.26243692885, 1189.12870543419,
    1239.75970033581, 1116.87300108878, 1077.66320677748, 1215.76397035317,
    1000.47341687018, 1138.8941620376, 977.591143870017, 1152.02164925958,
    969.223836265964, 1057.06097481807, 889.689253840196, 1165.81204298993,
    1216.08713830619, 1280.83371282797, 1208.4648457657, 1264.52605452378,
    1366.63902693168, 1191.6244929701, 1220.24603204681, 1323.03736956856,
    1090.30234711807, 1283.43746138696, 1157.43818887515, 1352.74337316786,
    1409.40984255262, 1476.26332329948, 1474.84261268216, 1532.03759868421,
    1377.6930108806, 1347.45715344042, 1312.27976150676, 1204.70172453162,
    1255.1321893586, 1240.03433129609, 1277.41968970551, 1340.36452560427,
    1248.51824770501, 1120.60434480781, 1177.3617450915, 1106.31929591054,
    1174.48548707999, 1100.56230467507, 781.054548693761, 658.656920660591,
    665.245385647703, 1436.1710767919, 1010.9713244342, 1085.95222996476,
    916.086429582412, 1005.67017841008, 1426.34551764972, 1378.64437398018,
    1375.90340420956, 1440.00216145498, 1178.04191915376, 1106.04723944121,
    774.579396923868, 821.079304972821, 856.59391417469, 862.385822263836,
    695.611212022243, 698.508436718344, 831.034863659039, 927.065683593056,
    800.244060459641, 861.81434371029, 1533.69752612507, 1549.59115191932,
    1507.27537337429, 1465.81955779263, 1347.78189536244, 1578.54212056382,
    1594.50099291575, 1058.37027757294, 1427.37717958346, 1431.27626895891,
    1362.62136841985, 1295.89537541058, 1249.52693245553, 1225.79184509631,
    1151.21572838259, 1156.7293492175, 1306.40107720521, 1142.26709470257,
    749.4102607034, 717.963831328974, 812.500978068188, 893.974083191674,
    648.662861395227, 602.74060711171, 712.480702845475, 740.512200550228,
    805.014126637161, 680.325930615691, 799.804801135588, 958.319354735665,
    874.162843365825, 889.76045155813, 921.352543327871, 859.634576000878,
    808.537585975561, 926.042441541467, 1004.14073021744, 1008.5058751011,
    1057.61375980799, 1491.69124999565, 1406.42002505818, 1411.39509682517,
    1482.35457167629, 1470.07127466057, 1434.5892538239, 755.133298363565,
    894.906504938554, 975.281264580714, 1483.06800654332, 1529.40491951512,
    1281.23999576587, 1311.28082687749, 1253.46889259152, 1103.89247077447,
    1105.06750735005, 1035.57918405457, 1199.16265439761, 1148.78514418515,
    1049.19702376991, 1158.60019241195, 1162.59311948116, 1095.65611656659,
    878.058374135844, 788.182906897712, 853.043884305736, 832.264563691969,
    741.587645145414, 617.746184145786, 1069.13588820083, 1052.05447961362,
    962.429499985256, 1028.66580320504, 1061.68257785881, 994.801286091748,
    924.628821793981, 909.406149524388, 991.153522095404, 1085.05250570571,
    969.157782621974, 936.087467324743, 1415.1359010846, 854.260557538138,
    925.653279962582 };

    std::vector<GridGeom::Edge> edges{
        {	142	,	184		}	,
        {	184	,	185		}	,
        {	142	,	185		}	,
        {	35	,	36		}	,
        {	36	,	185		}	,
        {	36	,	184		}	,
        {	37	,	184		}	,
        {	36	,	37		}	,
        {	38	,	228		}	,
        {	37	,	38		}	,
        {	37	,	228		}	,
        {	142	,	186		}	,
        {	142	,	231		}	,
        {	186	,	231		}	,
        {	37	,	186		}	,
        {	184	,	186		}	,
        {	186	,	228		}	,
        {	35	,	60		}	,
        {	36	,	60		}	,
        {	35	,	58		}	,
        {	58	,	60		}	,
        {	60	,	185		}	,
        {	60	,	143		}	,
        {	59	,	143		}	,
        {	59	,	60		}	,
        {	229	,	232		}	,
        {	59	,	229		}	,
        {	59	,	232		}	,
        {	34	,	61		}	,
        {	58	,	61		}	,
        {	34	,	58		}	,
        {	59	,	185		}	,
        {	142	,	229		}	,
        {	185	,	229		}	,
        {	229	,	231		}	,
        {	229	,	230		}	,
        {	230	,	231		}	,
        {	192	,	230		}	,
        {	192	,	229		}	,
        {	186	,	187		}	,
        {	187	,	231		}	,
        {	58	,	143		}	,
        {	34	,	35		}	,
        {	187	,	230		}	,
        {	39	,	148		}	,
        {	67	,	148		}	,
        {	39	,	67		}	,
        {	38	,	39		}	,
        {	38	,	148		}	,
        {	67	,	146		}	,
        {	39	,	146		}	,
        {	67	,	149		}	,
        {	48	,	67		}	,
        {	48	,	149		}	,
        {	148	,	228		}	,
        {	67	,	228		}	,
        {	146	,	149		}	,
        {	40	,	147		}	,
        {	41	,	147		}	,
        {	40	,	41		}	,
        {	40	,	146		}	,
        {	146	,	147		}	,
        {	66	,	146		}	,
        {	66	,	147		}	,
        {	66	,	149		}	,
        {	39	,	40		}	,
        {	147	,	155		}	,
        {	66	,	155		}	,
        {	81	,	149		}	,
        {	66	,	81		}	,
        {	48	,	248		}	,
        {	67	,	248		}	,
        {	186	,	248		}	,
        {	228	,	248		}	,
        {	187	,	244		}	,
        {	244	,	248		}	,
        {	187	,	248		}	,
        {	198	,	244		}	,
        {	187	,	198		}	,
        {	48	,	244		}	,
        {	48	,	81		}	,
        {	48	,	243		}	,
        {	81	,	243		}	,
        {	223	,	243		}	,
        {	224	,	243		}	,
        {	223	,	224		}	,
        {	48	,	224		}	,
        {	224	,	244		}	,
        {	81	,	223		}	,
        {	239	,	244		}	,
        {	224	,	239		}	,
        {	81	,	155		}	,
        {	198	,	230		}	,
        {	32	,	33		}	,
        {	32	,	61		}	,
        {	33	,	61		}	,
        {	31	,	32		}	,
        {	30	,	31		}	,
        {	32	,	62		}	,
        {	61	,	62		}	,
        {	47	,	61		}	,
        {	47	,	58		}	,
        {	47	,	143		}	,
        {	47	,	63		}	,
        {	63	,	143		}	,
        {	62	,	63		}	,
        {	62	,	188		}	,
        {	63	,	188		}	,
        {	47	,	62		}	,
        {	63	,	232		}	,
        {	143	,	232		}	,
        {	32	,	189		}	,
        {	62	,	189		}	,
        {	33	,	34		}	,
        {	144	,	188		}	,
        {	188	,	189		}	,
        {	144	,	189		}	,
        {	31	,	233		}	,
        {	189	,	233		}	,
        {	31	,	189		}	,
        {	30	,	233		}	,
        {	144	,	233		}	,
        {	144	,	190		}	,
        {	144	,	191		}	,
        {	190	,	191		}	,
        {	188	,	190		}	,
        {	191	,	200		}	,
        {	85	,	191		}	,
        {	85	,	200		}	,
        {	191	,	193		}	,
        {	144	,	193		}	,
        {	193	,	233		}	,
        {	191	,	194		}	,
        {	190	,	194		}	,
        {	30	,	193		}	,
        {	63	,	190		}	,
        {	190	,	232		}	,
        {	192	,	196		}	,
        {	196	,	230		}	,
        {	192	,	232		}	,
        {	196	,	198		}	,
        {	194	,	196		}	,
        {	192	,	194		}	,
        {	190	,	192		}	,
        {	195	,	196		}	,
        {	195	,	198		}	,
        {	227	,	239		}	,
        {	204	,	239		}	,
        {	204	,	227		}	,
        {	224	,	227		}	,
        {	195	,	204		}	,
        {	195	,	238		}	,
        {	204	,	238		}	,
        {	195	,	239		}	,
        {	198	,	239		}	,
        {	181	,	227		}	,
        {	181	,	204		}	,
        {	196	,	201		}	,
        {	195	,	201		}	,
        {	199	,	201		}	,
        {	196	,	199		}	,
        {	194	,	199		}	,
        {	194	,	200		}	,
        {	199	,	200		}	,
        {	197	,	200		}	,
        {	197	,	199		}	,
        {	197	,	247		}	,
        {	200	,	247		}	,
        {	197	,	201		}	,
        {	195	,	203		}	,
        {	203	,	238		}	,
        {	201	,	203		}	,
        {	203	,	234		}	,
        {	234	,	238		}	,
        {	202	,	203		}	,
        {	201	,	202		}	,
        {	202	,	234		}	,
        {	181	,	238		}	,
        {	197	,	245		}	,
        {	201	,	245		}	,
        {	139	,	238		}	,
        {	139	,	234		}	,
        {	137	,	234		}	,
        {	137	,	139		}	,
        {	180	,	227		}	,
        {	180	,	224		}	,
        {	41	,	42		}	,
        {	41	,	154		}	,
        {	42	,	154		}	,
        {	75	,	154		}	,
        {	77	,	154		}	,
        {	75	,	77		}	,
        {	42	,	75		}	,
        {	42	,	69		}	,
        {	42	,	49		}	,
        {	49	,	69		}	,
        {	43	,	69		}	,
        {	42	,	43		}	,
        {	74	,	75		}	,
        {	49	,	74		}	,
        {	49	,	75		}	,
        {	73	,	74		}	,
        {	73	,	75		}	,
        {	70	,	74		}	,
        {	49	,	70		}	,
        {	41	,	77		}	,
        {	43	,	44		}	,
        {	73	,	77		}	,
        {	65	,	77		}	,
        {	65	,	73		}	,
        {	77	,	147		}	,
        {	65	,	147		}	,
        {	65	,	81		}	,
        {	65	,	155		}	,
        {	80	,	81		}	,
        {	65	,	80		}	,
        {	79	,	80		}	,
        {	79	,	81		}	,
        {	73	,	76		}	,
        {	73	,	78		}	,
        {	76	,	78		}	,
        {	65	,	78		}	,
        {	80	,	83		}	,
        {	79	,	83		}	,
        {	78	,	83		}	,
        {	78	,	80		}	,
        {	74	,	76		}	,
        {	79	,	82		}	,
        {	82	,	83		}	,
        {	76	,	83		}	,
        {	1	,	44		}	,
        {	44	,	69		}	,
        {	1	,	69		}	,
        {	68	,	70		}	,
        {	2	,	70		}	,
        {	2	,	68		}	,
        {	68	,	69		}	,
        {	49	,	68		}	,
        {	1	,	68		}	,
        {	1	,	2		}	,
        {	3	,	71		}	,
        {	2	,	3		}	,
        {	2	,	71		}	,
        {	3	,	72		}	,
        {	71	,	72		}	,
        {	70	,	71		}	,
        {	3	,	4		}	,
        {	70	,	72		}	,
        {	72	,	150		}	,
        {	3	,	150		}	,
        {	72	,	83		}	,
        {	72	,	76		}	,
        {	151	,	210		}	,
        {	151	,	152		}	,
        {	152	,	210		}	,
        {	82	,	151		}	,
        {	72	,	151		}	,
        {	72	,	82		}	,
        {	150	,	151		}	,
        {	4	,	150		}	,
        {	150	,	210		}	,
        {	4	,	210		}	,
        {	5	,	210		}	,
        {	4	,	5		}	,
        {	5	,	153		}	,
        {	153	,	210		}	,
        {	5	,	6		}	,
        {	82	,	152		}	,
        {	70	,	76		}	,
        {	79	,	182		}	,
        {	82	,	182		}	,
        {	222	,	223		}	,
        {	180	,	223		}	,
        {	180	,	222		}	,
        {	79	,	223		}	,
        {	226	,	227		}	,
        {	180	,	226		}	,
        {	179	,	180		}	,
        {	179	,	226		}	,
        {	79	,	178		}	,
        {	178	,	182		}	,
        {	79	,	222		}	,
        {	178	,	222		}	,
        {	64	,	182		}	,
        {	64	,	178		}	,
        {	178	,	179		}	,
        {	179	,	222		}	,
        {	177	,	178		}	,
        {	177	,	179		}	,
        {	181	,	226		}	,
        {	179	,	181		}	,
        {	152	,	182		}	,
        {	136	,	181		}	,
        {	181	,	225		}	,
        {	136	,	225		}	,
        {	139	,	181		}	,
        {	139	,	225		}	,
        {	138	,	225		}	,
        {	131	,	138		}	,
        {	131	,	225		}	,
        {	138	,	139		}	,
        {	131	,	136		}	,
        {	136	,	177		}	,
        {	135	,	136		}	,
        {	135	,	177		}	,
        {	136	,	179		}	,
        {	134	,	135		}	,
        {	129	,	134		}	,
        {	129	,	135		}	,
        {	134	,	136		}	,
        {	131	,	134		}	,
        {	64	,	177		}	,
        {	131	,	133		}	,
        {	133	,	134		}	,
        {	64	,	135		}	,
        {	145	,	152		}	,
        {	145	,	153		}	,
        {	152	,	153		}	,
        {	64	,	152		}	,
        {	64	,	145		}	,
        {	64	,	206		}	,
        {	145	,	206		}	,
        {	135	,	206		}	,
        {	206	,	209		}	,
        {	145	,	209		}	,
        {	153	,	205		}	,
        {	5	,	205		}	,
        {	145	,	205		}	,
        {	6	,	209		}	,
        {	7	,	209		}	,
        {	6	,	7		}	,
        {	6	,	205		}	,
        {	205	,	209		}	,
        {	206	,	208		}	,
        {	208	,	209		}	,
        {	169	,	207		}	,
        {	169	,	208		}	,
        {	207	,	208		}	,
        {	207	,	246		}	,
        {	169	,	246		}	,
        {	129	,	246		}	,
        {	129	,	207		}	,
        {	130	,	134		}	,
        {	129	,	130		}	,
        {	135	,	207		}	,
        {	206	,	207		}	,
        {	214	,	215		}	,
        {	169	,	215		}	,
        {	169	,	214		}	,
        {	7	,	208		}	,
        {	208	,	215		}	,
        {	7	,	215		}	,
        {	8	,	215		}	,
        {	7	,	8		}	,
        {	8	,	214		}	,
        {	214	,	246		}	,
        {	128	,	129		}	,
        {	128	,	246		}	,
        {	8	,	9		}	,
        {	85	,	193		}	,
        {	29	,	30		}	,
        {	29	,	160		}	,
        {	30	,	160		}	,
        {	161	,	211		}	,
        {	160	,	211		}	,
        {	160	,	161		}	,
        {	85	,	160		}	,
        {	160	,	193		}	,
        {	245	,	247		}	,
        {	165	,	247		}	,
        {	85	,	247		}	,
        {	85	,	165		}	,
        {	85	,	211		}	,
        {	165	,	211		}	,
        {	86	,	161		}	,
        {	86	,	211		}	,
        {	86	,	165		}	,
        {	29	,	161		}	,
        {	28	,	161		}	,
        {	28	,	87		}	,
        {	87	,	161		}	,
        {	28	,	29		}	,
        {	156	,	162		}	,
        {	87	,	162		}	,
        {	87	,	156		}	,
        {	27	,	87		}	,
        {	27	,	28		}	,
        {	87	,	164		}	,
        {	161	,	164		}	,
        {	86	,	164		}	,
        {	164	,	165		}	,
        {	50	,	164		}	,
        {	162	,	164		}	,
        {	50	,	162		}	,
        {	162	,	240		}	,
        {	50	,	240		}	,
        {	84	,	156		}	,
        {	84	,	162		}	,
        {	163	,	212		}	,
        {	162	,	163		}	,
        {	162	,	212		}	,
        {	27	,	156		}	,
        {	50	,	165		}	,
        {	236	,	245		}	,
        {	202	,	236		}	,
        {	202	,	245		}	,
        {	165	,	245		}	,
        {	165	,	236		}	,
        {	165	,	241		}	,
        {	236	,	241		}	,
        {	50	,	241		}	,
        {	236	,	237		}	,
        {	213	,	236		}	,
        {	213	,	237		}	,
        {	137	,	235		}	,
        {	234	,	235		}	,
        {	235	,	236		}	,
        {	202	,	235		}	,
        {	120	,	141		}	,
        {	141	,	237		}	,
        {	120	,	237		}	,
        {	141	,	235		}	,
        {	235	,	237		}	,
        {	137	,	141		}	,
        {	213	,	241		}	,
        {	213	,	221		}	,
        {	213	,	240		}	,
        {	221	,	240		}	,
        {	240	,	241		}	,
        {	163	,	242		}	,
        {	240	,	242		}	,
        {	163	,	240		}	,
        {	106	,	242		}	,
        {	106	,	163		}	,
        {	221	,	242		}	,
        {	120	,	221		}	,
        {	220	,	221		}	,
        {	120	,	220		}	,
        {	120	,	213		}	,
        {	107	,	220		}	,
        {	173	,	220		}	,
        {	107	,	173		}	,
        {	173	,	221		}	,
        {	173	,	242		}	,
        {	120	,	183		}	,
        {	141	,	183		}	,
        {	106	,	173		}	,
        {	120	,	122		}	,
        {	122	,	220		}	,
        {	84	,	212		}	,
        {	26	,	156		}	,
        {	26	,	157		}	,
        {	156	,	157		}	,
        {	26	,	27		}	,
        {	157	,	158		}	,
        {	25	,	157		}	,
        {	25	,	158		}	,
        {	84	,	110		}	,
        {	84	,	157		}	,
        {	110	,	157		}	,
        {	106	,	212		}	,
        {	110	,	212		}	,
        {	106	,	110		}	,
        {	94	,	158		}	,
        {	94	,	110		}	,
        {	110	,	158		}	,
        {	94	,	159		}	,
        {	158	,	159		}	,
        {	108	,	110		}	,
        {	94	,	108		}	,
        {	25	,	26		}	,
        {	25	,	159		}	,
        {	24	,	25		}	,
        {	24	,	159		}	,
        {	96	,	159		}	,
        {	24	,	96		}	,
        {	23	,	24		}	,
        {	24	,	54		}	,
        {	23	,	54		}	,
        {	46	,	54		}	,
        {	53	,	54		}	,
        {	46	,	53		}	,
        {	23	,	53		}	,
        {	22	,	53		}	,
        {	22	,	45		}	,
        {	45	,	53		}	,
        {	21	,	22		}	,
        {	22	,	23		}	,
        {	54	,	96		}	,
        {	94	,	96		}	,
        {	21	,	45		}	,
        {	106	,	109		}	,
        {	106	,	108		}	,
        {	108	,	109		}	,
        {	96	,	104		}	,
        {	55	,	96		}	,
        {	55	,	104		}	,
        {	104	,	108		}	,
        {	94	,	104		}	,
        {	104	,	109		}	,
        {	105	,	219		}	,
        {	109	,	219		}	,
        {	105	,	109		}	,
        {	173	,	219		}	,
        {	109	,	173		}	,
        {	101	,	102		}	,
        {	55	,	102		}	,
        {	55	,	101		}	,
        {	102	,	109		}	,
        {	102	,	104		}	,
        {	102	,	105		}	,
        {	95	,	96		}	,
        {	55	,	95		}	,
        {	57	,	95		}	,
        {	46	,	95		}	,
        {	46	,	57		}	,
        {	54	,	95		}	,
        {	46	,	56		}	,
        {	45	,	46		}	,
        {	45	,	56		}	,
        {	20	,	45		}	,
        {	20	,	21		}	,
        {	19	,	45		}	,
        {	19	,	20		}	,
        {	55	,	57		}	,
        {	57	,	101		}	,
        {	57	,	97		}	,
        {	56	,	57		}	,
        {	56	,	97		}	,
        {	19	,	56		}	,
        {	19	,	97		}	,
        {	18	,	19		}	,
        {	18	,	97		}	,
        {	97	,	98		}	,
        {	57	,	98		}	,
        {	107	,	219		}	,
        {	137	,	138		}	,
        {	137	,	140		}	,
        {	138	,	140		}	,
        {	140	,	141		}	,
        {	133	,	140		}	,
        {	133	,	138		}	,
        {	128	,	130		}	,
        {	130	,	133		}	,
        {	132	,	133		}	,
        {	130	,	132		}	,
        {	119	,	130		}	,
        {	119	,	132		}	,
        {	119	,	128		}	,
        {	132	,	140		}	,
        {	89	,	128		}	,
        {	89	,	246		}	,
        {	117	,	140		}	,
        {	117	,	132		}	,
        {	117	,	183		}	,
        {	118	,	183		}	,
        {	117	,	118		}	,
        {	117	,	141		}	,
        {	113	,	118		}	,
        {	114	,	118		}	,
        {	113	,	114		}	,
        {	118	,	122		}	,
        {	122	,	183		}	,
        {	118	,	121		}	,
        {	117	,	121		}	,
        {	119	,	121		}	,
        {	121	,	123		}	,
        {	119	,	123		}	,
        {	121	,	132		}	,
        {	51	,	123		}	,
        {	116	,	123		}	,
        {	51	,	116		}	,
        {	113	,	123		}	,
        {	113	,	121		}	,
        {	119	,	124		}	,
        {	124	,	128		}	,
        {	114	,	122		}	,
        {	123	,	124		}	,
        {	127	,	214		}	,
        {	89	,	214		}	,
        {	89	,	127		}	,
        {	92	,	126		}	,
        {	125	,	126		}	,
        {	92	,	125		}	,
        {	124	,	126		}	,
        {	89	,	126		}	,
        {	89	,	124		}	,
        {	126	,	127		}	,
        {	8	,	127		}	,
        {	9	,	127		}	,
        {	127	,	172		}	,
        {	9	,	172		}	,
        {	92	,	172		}	,
        {	92	,	171		}	,
        {	171	,	172		}	,
        {	92	,	127		}	,
        {	10	,	172		}	,
        {	9	,	10		}	,
        {	124	,	125		}	,
        {	52	,	125		}	,
        {	52	,	171		}	,
        {	125	,	171		}	,
        {	125	,	168		}	,
        {	52	,	168		}	,
        {	51	,	168		}	,
        {	51	,	125		}	,
        {	88	,	168		}	,
        {	51	,	88		}	,
        {	51	,	124		}	,
        {	11	,	52		}	,
        {	52	,	167		}	,
        {	11	,	167		}	,
        {	10	,	171		}	,
        {	91	,	171		}	,
        {	10	,	91		}	,
        {	11	,	91		}	,
        {	10	,	11		}	,
        {	52	,	91		}	,
        {	167	,	168		}	,
        {	113	,	116		}	,
        {	111	,	219		}	,
        {	107	,	111		}	,
        {	107	,	122		}	,
        {	105	,	111		}	,
        {	103	,	105		}	,
        {	103	,	111		}	,
        {	111	,	112		}	,
        {	107	,	112		}	,
        {	115	,	116		}	,
        {	113	,	115		}	,
        {	114	,	115		}	,
        {	112	,	114		}	,
        {	112	,	115		}	,
        {	90	,	112		}	,
        {	90	,	115		}	,
        {	107	,	114		}	,
        {	103	,	112		}	,
        {	90	,	103		}	,
        {	93	,	105		}	,
        {	93	,	103		}	,
        {	93	,	99		}	,
        {	93	,	101		}	,
        {	99	,	101		}	,
        {	93	,	102		}	,
        {	98	,	101		}	,
        {	98	,	99		}	,
        {	18	,	98		}	,
        {	17	,	98		}	,
        {	17	,	18		}	,
        {	217	,	218		}	,
        {	100	,	218		}	,
        {	100	,	217		}	,
        {	93	,	218		}	,
        {	93	,	100		}	,
        {	99	,	100		}	,
        {	16	,	99		}	,
        {	16	,	100		}	,
        {	17	,	99		}	,
        {	16	,	17		}	,
        {	103	,	218		}	,
        {	103	,	216		}	,
        {	216	,	218		}	,
        {	88	,	175		}	,
        {	51	,	175		}	,
        {	174	,	175		}	,
        {	116	,	175		}	,
        {	116	,	174		}	,
        {	90	,	176		}	,
        {	116	,	176		}	,
        {	90	,	116		}	,
        {	90	,	170		}	,
        {	170	,	176		}	,
        {	88	,	174		}	,
        {	88	,	167		}	,
        {	88	,	166		}	,
        {	166	,	167		}	,
        {	12	,	166		}	,
        {	13	,	166		}	,
        {	12	,	13		}	,
        {	12	,	167		}	,
        {	13	,	88		}	,
        {	11	,	12		}	,
        {	174	,	176		}	,
        {	170	,	174		}	,
        {	170	,	216		}	,
        {	216	,	217		}	,
        {	170	,	217		}	,
        {	15	,	217		}	,
        {	15	,	170		}	,
        {	16	,	217		}	,
        {	15	,	16		}	,
        {	13	,	174		}	,
        {	14	,	174		}	,
        {	13	,	14		}	,
        {	14	,	170		}	,
        {	14	,	15		}	,
        {	90	,	216		} };

    std::vector<GridGeom::Point> nodes(xCoordinates.size());

    for (int i = 0; i < nodes.size(); i++)
    {
        nodes[i].x = xCoordinates[i];
        nodes[i].y = yCoordinates[i];
    }

    for (int i = 0; i < edges.size(); i++)
    {
        edges[i].first -= 1;
        edges[i].second -= 1;
    }

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 0;
    int projectToLandBoundaryOption = 0;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    GridGeomApi::GeometryListNative geometryListNativePolygon;
    GridGeomApi::GeometryListNative geometryListNativeLandBoundaries;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.5;

    // now build node-edge mapping
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);
    GridGeom::Orthogonalization orthogonalization;

    std::vector<GridGeom::Point> polygon;
    std::vector<GridGeom::Point> landBoundary;

    orthogonalization.Set(mesh,
        isTriangulationRequired,
        isAccountingForLandBoundariesRequired,
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);

    orthogonalization.Iterate(mesh);

    constexpr double tolerance = 1.5;

    // check the first 10 points
    ASSERT_NEAR(68.771705432835475, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(169.49338272334273, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(262.80128484924921, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(361.60010033352023, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(468.13991812406925, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(549.89461192844624, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(653.02704974527421, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(747.81537706979441, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(853.40641427112951, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(938.69752431820143, mesh.m_nodes[9].x, tolerance);

    ASSERT_NEAR(1399.7751472360221, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(1426.5945287630802, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(1451.4398281457179, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(1477.7472050498141, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(1506.1157955857589, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(1527.8847968946166, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(1555.3460969050145, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(1580.5855923464549, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(1608.7015489976982, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(1631.412199601948, mesh.m_nodes[9].y, tolerance);

}

TEST(OrthogonalizationTests, TestOrthogonalizationFourQuads)
{

    const int n = 3; //x
    const int m = 3; //y

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<GridGeom::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            indexesValues[i][j] = i + j * n;
            nodes[nodeIndex] = { (double)i, (double)j };
            nodeIndex++;
        }
    }

    std::vector<GridGeom::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = { indexesValues[i][j], indexesValues[i + 1][j] };
            edgeIndex++;
        }
    }

    for (int j = 0; j < m - 1; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            edges[edgeIndex] = { indexesValues[i][j + 1], indexesValues[i][j] };
            edgeIndex++;
        }
    }

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 0;
    int projectToLandBoundaryOption = 0;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;

    // now build node-edge mapping
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    std::vector<GridGeom::Point> polygon;
    std::vector<GridGeom::Point> landBoundary;

    GridGeom::Orthogonalization orthogonalization;
    orthogonalization.Set(mesh,
        isTriangulationRequired,
        isAccountingForLandBoundariesRequired,
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);
}

TEST(OrthogonalizationTests, OrthogonalizeAndSnapToLandBoundaries)
{
    // Prepare
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 322.252624511719,454.880187988281 });
    nodes.push_back(GridGeom::Point{ 227.002044677734,360.379241943359 });
    nodes.push_back(GridGeom::Point{ 259.252227783203,241.878051757813 });
    nodes.push_back(GridGeom::Point{ 428.003295898438,210.377746582031 });
    nodes.push_back(GridGeom::Point{ 536.003967285156,310.878753662109 });
    nodes.push_back(GridGeom::Point{ 503.753784179688,432.379974365234 });
    nodes.push_back(GridGeom::Point{ 350.752807617188,458.630249023438 });
    nodes.push_back(GridGeom::Point{ 343.15053976393,406.232256102912 });
    nodes.push_back(GridGeom::Point{ 310.300984548069,319.41005739802 });
    nodes.push_back(GridGeom::Point{ 423.569603308318,326.17986967523 });

    std::vector<GridGeom::Edge> edges;
    edges.push_back({ 3, 9 });
    edges.push_back({ 9, 2 });
    edges.push_back({ 2, 3 });
    edges.push_back({ 3, 4 });
    edges.push_back({ 4, 9 });
    edges.push_back({ 2, 8 });
    edges.push_back({ 8, 1 });
    edges.push_back({ 1, 2 });
    edges.push_back({ 9, 8 });
    edges.push_back({ 8, 7 });
    edges.push_back({ 7, 1 });
    edges.push_back({ 9, 10 });
    edges.push_back({ 10, 8 });
    edges.push_back({ 4, 5 });
    edges.push_back({ 5, 10 });
    edges.push_back({ 10, 4 });
    edges.push_back({ 8, 6 });
    edges.push_back({ 6, 7 });
    edges.push_back({ 10, 6 });
    edges.push_back({ 5, 6 });

    for (int i = 0; i < edges.size(); i++)
    {
        edges[i].first -= 1;
        edges[i].second -= 1;
    }

    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    // the land boundary to use
    std::vector<GridGeom::Point> landBoundary
    {
        { 235.561218, 290.571899 },
        { 265.953522, 436.515747 },
        { 429.349854, 450.959656 },
        { 535.271545, 386.262909 },
        { GridGeom::doubleMissingValue, GridGeom::doubleMissingValue },
        { 246.995941, 262.285858 },
        { 351.112183, 237.309906 },
        { 443.191895, 262.285858 },
        { 553.627319, 327.283539 },
    };

    // no enclosing polygon
    std::vector<GridGeom::Point> polygon;

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 1;

    // snap to land boundaries
    int projectToLandBoundaryOption = 2;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;

    GridGeom::Orthogonalization orthogonalization;
    orthogonalization.Set(mesh,
        isTriangulationRequired,
        isAccountingForLandBoundariesRequired,
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);

    orthogonalization.Iterate(mesh);

    // check the values
    constexpr double tolerance = 0.15;
    ASSERT_NEAR(313.081472564480, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(253.641466857330, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(254.777224294204, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(443.191895000000, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(535.240231516760, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(480.436129612752, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(345.948240805397, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(342.668434889472, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(318.414413615199, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(424.616311031376, mesh.m_nodes[9].x, tolerance);

    ASSERT_NEAR(440.681763586650, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(377.393256506700, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(260.419242817573, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(262.285858000000, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(316.461666783032, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(419.756265860671, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(443.587120174434, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(402.913858250569, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(336.831643075189, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(340.875100904741, mesh.m_nodes[9].y, tolerance);
}