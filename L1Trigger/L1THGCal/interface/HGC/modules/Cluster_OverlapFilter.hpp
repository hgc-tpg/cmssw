// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------
#pragma once

#include <array>
#include <cstdlib>

#include "L1Trigger/L1THGCal/interface/HGC/HgcTypes.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/generated/PkgLutValues.hpp"

namespace HGC
{

  void Cluster_OverlapFilter_Flags1( const std::array< std::array< Cluster , 44 > ,  41 >& ClusterIn ,
                                           std::array< std::array< DebugBoolean ,  44 > , 41 >& CentreSel ,
                                           std::array< std::array< DebugBoolean ,  44 > , 41 >& UpGtDn ,
                                           std::array< std::array< DebugBoolean ,  44 > , 41 >& UpGtDnLt ,
                                           std::array< std::array< DebugBoolean ,  44 > , 41 >& UpGtDnRt ){

    for ( int lCol(0); lCol!=44 ; ++lCol ) { // col
      for ( int lRow(0); lRow!=41 ; ++lRow ) { // row

        // Add the debugging-fields to the output cell
        const Cluster& Up = ClusterIn.at( lRow ).at( lCol );

        {
          auto& x = CentreSel.at( lRow ).at( lCol );
          x.A = Up;
          x.B = Cluster();
          x.Value = Up.Field0.E > 0;
        }

        {
          auto& x = UpGtDn.at( lRow ).at( lCol );
          x.A = Up;
          if( lRow > 1 ) {
            const auto& Dn = ClusterIn.at( lRow-2 ).at( lCol );
            x.B = Dn;
            x.Value = ( Up.Field0.E > Dn.Field0.E ) || !( Up.DataValid && Dn.DataValid );
          } else {
            x.B = Cluster();
            x.Value = true;
          }
        }

        {
          auto& x = UpGtDnLt.at( lRow ).at( lCol );
          x.A = Up;
          if( lRow > 0 and lCol > 0 ) {
            const auto& DnLt = ClusterIn.at( lRow-1 ).at( lCol-1 );
            x.B = DnLt;
            x.Value = ( Up.Field0.E > DnLt.Field0.E ) || !( Up.DataValid && DnLt.DataValid );
          } else {
            x.B = Cluster();
            x.Value = true;
          }
        }

        {
          auto& x = UpGtDnRt.at( lRow ).at( lCol );
          x.A = Up;
          if( lRow > 0 and lCol < 43 ) {
            const auto& DnRt = ClusterIn.at( lRow-1 ).at( lCol+1 );
            x.B = DnRt;
            x.Value = ( Up.Field0.E > DnRt.Field0.E ) || !( Up.DataValid && DnRt.DataValid );
          } else {
            x.B = Cluster();
            x.Value = true;
          }
        }

      }
    }

  }


  void Cluster_OverlapFilter( const std::array< std::array< Cluster       , 44 > , 41 >& ClusterIn ,
                                    std::array< std::array< DebugBoolean2 , 44 > , 41 >& DebugOut ,
                                    std::array< std::array< Cluster       , 44 > , 41 >& ClusterOut )
  {

    ClusterOut = std::array< std::array< Cluster , 44 > ,  41 >(); // Start afresh
    DebugOut   = std::array< std::array< DebugBoolean2 , 44 > ,  41 >(); // Start afresh

    for ( int lCol(0); lCol!=44 ; ++lCol ) { // col
      for ( int lRow(0); lRow!=41 ; ++lRow ) { // row
        // Add the debugging-fields to the output cell
        const Cluster& lIn = ClusterIn.at( lRow ).at( lCol );

        Cluster& lOut        = ClusterOut.at( lRow ).at( lCol );
        DebugBoolean2& lDebug = DebugOut.at( lRow ).at( lCol );

        lOut = lIn;
        lOut.ColumnSet = ( (lCol/4) - 5 ) & 0xF;

        // centre
        lOut.DataValid &= ( lOut.Field0.E > 0 );

        lDebug.Ce = lIn;

        // up
        if( lRow < 38 ){
          lOut.DataValid &= ( lIn.Field0.E >= ClusterIn.at( lRow+2 ).at( lCol ).Field0.E );
          lDebug.Up = ClusterIn.at( lRow+2 ).at( lCol );
        }

        if( lCol < 43 ){
          // up right
          if( lRow < 39 ){
            lOut.DataValid &= ( lIn.Field0.E >= ClusterIn.at( lRow+1 ).at( lCol+1 ).Field0.E );
            lDebug.UpRt = ClusterIn.at( lRow+1 ).at( lCol+1 );
          }
          // down right
          if( lRow >  0 ){
            lOut.DataValid &= ( lIn.Field0.E >= ClusterIn.at( lRow-1 ).at( lCol+1 ).Field0.E );
            lDebug.DnRt = ClusterIn.at( lRow-1 ).at( lCol+1 );
          }
        }

        // down
        if( lRow > 1 ){
          lOut.DataValid &= ( lIn.Field0.E >  ClusterIn.at( lRow-2 ).at( lCol ).Field0.E );
          lDebug.Dn = ClusterIn.at( lRow-2 ).at( lCol );
        }

        if( lCol > 0 ){
          // down left
          if( lRow >  0 ){
            lOut.DataValid &= ( lIn.Field0.E >  ClusterIn.at( lRow-1 ).at( lCol-1 ).Field0.E );
            lDebug.DnLt = ClusterIn.at( lRow-1 ).at( lCol-1 );
          }
          // up left
          if( lRow < 39 ){
            lOut.DataValid &= ( lIn.Field0.E >  ClusterIn.at( lRow+1 ).at( lCol-1 ).Field0.E );
            lDebug.UpLt = ClusterIn.at( lRow+1 ).at( lCol-1 );
          }
        }


        lDebug.Value = lOut.DataValid;
      }
    }


  }

}
