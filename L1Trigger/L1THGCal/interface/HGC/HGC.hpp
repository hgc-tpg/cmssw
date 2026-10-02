

#include "L1Trigger/L1THGCal/interface/HGC/modules/Links_UnpackInputs.hpp"

#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_Routing.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_Decoders.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_Accumulator.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_ColumnAdder.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_HexagonSum.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_OverlapFilter.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_TriangleFilter.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_Funnel.hpp"
// #include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_Buffer.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_Buffer2.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_Properties.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Cluster_PackLinks.hpp"

#include "L1Trigger/L1THGCal/interface/HGC/modules/Tower_Unpack.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Tower_Sum.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Tower_Format.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/modules/Tower_Buffer.hpp"

#include "L1Trigger/L1THGCal/interface/HGC/modules/Links_MuxOutput.hpp"

#include <memory>

namespace HGC
{

  // ===================================================================================================
  __attribute__((flatten)) // Inline all content
  void Clusters_Step1( const std::array< std::array< LinkTriggerCell , 192 > , 162 >& LinkCellsIn ,
                             std::array< std::array< Cluster         , 154 > , 162 >& ProtoClustersOut )
  {
    // Intermediate signals for chaining algorithm modules
    std::array< std::array< LinkTriggerCell , 289 > , 162 > RoutedCellsInt; // Intermediate
    std::array< std::array< LinkTriggerCell , 154 > , 162 > RoutedCellsOut;

    // Emulate the algorithm modules
    Cluster_Routing( LinkCellsIn , RoutedCellsInt , RoutedCellsOut );
    Cluster_Decoders( RoutedCellsOut , ProtoClustersOut );
  }
  // ===================================================================================================

  // ===================================================================================================
  __attribute__((flatten)) // Inline all content
  void Clusters_Step2( const std::array< std::array< Cluster , 154 > , 162 >& ProtoClustersIn,
                             std::array< ClusterProperty             , 162 >& ClusterPropertiesOut )
  {
    // These signals occupy roughly 7.5 MiB.  Keep them off the CMSSW worker
    // thread stack, which is normally limited to 8 MiB and also contains the
    // input and output arrays owned by the caller.
    struct Workspace {
      std::array< std::array< Cluster , 154 > , 162 > AccumulatedClusters;
      std::array< std::array< Cluster ,  44 > ,  41 > Triangles1 , Hexagons1 , Filtered1;
      std::array< std::array< Cluster ,  44 > ,  41 > Triangles2 , Hexagons2 , Filtered2;
      std::array< std::array< Cluster ,  44 > ,  41 > Triangles3;
      std::array< std::array< DebugBoolean2 , 44 > , 41 > Debug1 , Debug2;
      std::array< Cluster ,  162 > Funneled1 , Funneled2 , Buffered;
    };
    auto workspace = std::make_unique<Workspace>();

    // Emulate the algorithm modules
    Cluster_Accumulator( ProtoClustersIn , workspace->AccumulatedClusters );
    Cluster_ColumnAdder( workspace->AccumulatedClusters , workspace->Triangles1 );
    Cluster_HexagonSum( workspace->Triangles1 , workspace->Hexagons1 );
    Cluster_OverlapFilter( workspace->Hexagons1 , workspace->Debug1 , workspace->Filtered1 );
    Cluster_TriangleFilter( workspace->Hexagons1 , workspace->Triangles1 , workspace->Triangles2 );
    Cluster_Funnel( workspace->Hexagons1 , workspace->Funneled1 );
    Cluster_HexagonSum( workspace->Triangles2 , workspace->Hexagons2 );
    Cluster_OverlapFilter( workspace->Hexagons2 , workspace->Debug2 , workspace->Filtered2 );
    Cluster_TriangleFilter( workspace->Hexagons2 , workspace->Triangles2 , workspace->Triangles3 );
    Cluster_Funnel( workspace->Hexagons2 , workspace->Funneled2 );
    Cluster_Buffer2( workspace->Funneled1 , workspace->Funneled2 , workspace->Buffered );
    Cluster_Properties( workspace->Buffered , ClusterPropertiesOut );
  }
  // ===================================================================================================

  // ===================================================================================================
  __attribute__((flatten)) // Inline all content
  void Clusters( const std::array< std::array< LinkTriggerCell , 192 > , 162 >& LinkCellsIn ,
                       std::array< ClusterProperty                     , 162 >& ClusterPropertiesOut )
  {

    // Intermediate signals for chaining algorithm modules
    std::array< std::array< Cluster , 154 > , 162 > ProtoClusters;

    // Emulate the algorithm modules
    Clusters_Step1( LinkCellsIn     , ProtoClusters     );
    Clusters_Step2( ProtoClusters , ClusterPropertiesOut );
  }
  // ===================================================================================================

  // ===================================================================================================
  __attribute__((flatten)) // Inline all content
  void Clusters( const std::array< std::array< LinkTriggerCell , 192 > , 162 >& LinkCellsIn ,
                       std::array< std::array< lword     ,  4 > , 162 >& ClusterLinksOut )
  {

    // Intermediate signals for chaining algorithm modules
    std::array< ClusterProperty             , 162 > ClusterProperties;

    // Emulate the algorithm modules
    Clusters( LinkCellsIn     , ClusterProperties );
    Cluster_PackLinks( ClusterProperties , ClusterLinksOut );
  }
  // ===================================================================================================


  // ===================================================================================================
  __attribute__((flatten)) // Inline all content
  void Towers( const std::array< std::array< LinkTower , 96 > , 162 >& TowersIn ,
                     std::array< std::array< FormattedTower , 24 > , 20 >& FormattedTowersOut )
  {
    // Intermediate signals for chaining algorithm modules
    std::array< std::array< Tower          , 96 > , 162 > UnpackedTowers;
    std::array< std::array< Tower          ,  5 > , 162 > SummedTowers;
    std::array< std::array< FormattedTower ,  5 > , 162 > FormattedTowers;
    std::array< std::array< PackedTower    ,  4 > , 162 > TowerDebug;
    std::array< std::array< lword          ,  4 > , 162 > TowerLinksOut;

    // Emulate the algorithm modules
    Tower_Unpack( TowersIn , UnpackedTowers );
    Tower_Sum( UnpackedTowers , SummedTowers );
    Tower_Format( SummedTowers , FormattedTowers );
    Tower_Buffer( FormattedTowers , FormattedTowersOut , TowerLinksOut );
  }
  // ===================================================================================================

  // ===================================================================================================
  __attribute__((flatten)) // Inline all content
  void Towers( const std::array< std::array< LinkTower , 96 > , 162 >& TowersIn ,
                     std::array< std::array< lword     ,  4 > , 162 >& TowerLinksOut )
  {
    // Intermediate signals for chaining algorithm modules
    std::array< std::array< Tower          , 96 > , 162 > UnpackedTowers;
    std::array< std::array< Tower          ,  5 > , 162 > SummedTowers;
    std::array< std::array< FormattedTower ,  5 > , 162 > FormattedTowers;
    std::array< std::array< PackedTower    ,  4 > , 162 > TowerDebug;
    std::array< std::array< FormattedTower , 24 > ,  20 > FormattedTowersOut;

    // Emulate the algorithm modules
    Tower_Unpack( TowersIn , UnpackedTowers );
    Tower_Sum( UnpackedTowers , SummedTowers );
    Tower_Format( SummedTowers , FormattedTowers );
    Tower_Buffer( FormattedTowers , FormattedTowersOut , TowerLinksOut );
  }
  // ===================================================================================================


  // ===================================================================================================
  __attribute__((flatten)) // Inline all content
  __attribute__((warning("Function 'Stage2' is an example only, lacking the necessary configuration to make it meaningful")))
  void Stage2( const std::array< std::array< lword , 96 > , 162 >& LinksIn ,
                     std::array< std::array< lword ,  4 > , 162 >& LinksOut )
  {

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!
    // !!! EXAMPLE ONLY
    // !!! This will not work correctly because the latency of the tower buffer must be tuned
    // !!! to match the latency of the clusters. Currently there is no way to do so.
    // !!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // Intermediate signals for chaining algorithm modules
    std::array< std::array< LinkTriggerCell , 96*2 > , 162 > LinkCells;
    std::array< std::array< LinkTower       , 96   > , 162 > LinkTowers;
    std::array< std::array< LinkFlags       , 96   > , 162 > LinkFlags;
    std::array< std::array< lword           ,    4 > , 162 > TowerLinks;
    std::array< std::array< lword           ,    4 > , 162 > ClusterLinks;

    // Emulate the algorithm modules
    Links_UnpackInput( LinksIn , LinkCells , LinkTowers , LinkFlags );
    Towers( LinkTowers , TowerLinks );
    Clusters( LinkCells , ClusterLinks );
    Links_MuxOutput( TowerLinks , ClusterLinks , LinksOut );

  }
  // ===================================================================================================



}
