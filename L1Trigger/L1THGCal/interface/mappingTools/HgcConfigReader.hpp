#pragma once

#include <tuple>
#include <unordered_map>
#include <cstdint>
#include <string>
#include <vector>

namespace l1thgcmapping {

  // Commented out rather than removed as will be required for firmware mapping
  // However may end up in another class/file
  // struct tuplehash {
  // private:

  //   template< typename Tuple, size_t... I >
  //   std::size_t combine_hash( const Tuple& t, std::index_sequence<I...> ) const
  //   {
  //     constexpr std::size_t size = sizeof...(I);
  //     auto seed = combine_hash( t , std::make_index_sequence<size-1>{} );
  //     auto current = std::get<size-1>( t );
  //     std::hash< decltype( current ) > hasher;
  //     return seed ^ ( hasher( current ) + 0x9e3779b9 + (seed<<6) + (seed>>2) ); // Adapted from Boost Hash_Combine
  //   }

  //   template< typename Tuple >
  //   std::size_t combine_hash( const Tuple& t, std::index_sequence<> ) const
  //   {
  //     return 0;
  //   }

  // public:    
  //   template< typename... Args > std::size_t operator()( const std::tuple< Args... > &x ) const { return combine_hash( x, std::make_index_sequence<sizeof...(Args)>{}); }
    
  // };



  struct Maps
  {
    typedef uint32_t tID;
    typedef std::tuple< tID , tID > tID2;   
    typedef std::tuple< tID , tID , tID > tID3;    
    typedef std::unordered_map < tID , tID > tId_to_tId;
    typedef std::unordered_map < tID , std::vector<tID> > tId_to_tIds;
    typedef std::unordered_multimap < tID , tID > multi_tId_to_tId;
    // typedef std::unordered_map < tID2 , tID3 , tuplehash > tId2_to_tId3;  
    // typedef std::unordered_map < tID3 , tID2 , tuplehash > tId3_to_tId2;  
    
    tId_to_tId       module_to_motherboard;
    multi_tId_to_tId motherboard_to_module;

    tId_to_tId       motherboard_to_region;
    multi_tId_to_tId region_to_motherboard;  

    tId_to_tId       region_to_stage1;
    multi_tId_to_tId stage1_to_region;  

    // tId3_to_tId2     moduleColIndex_to_s1ChannelFrame;
    // tId2_to_tId3     s1ChannelFrame_to_moduleColIndex; 
    
    tId_to_tId       belink_to_reflink;

    tId_to_tId       s1out_to_belink;
    tId_to_tId       belink_to_s1out;

    tId_to_tId       belink_to_s2in;
    tId_to_tId       s2in_to_belink; 
    
    tId_to_tIds      motherboard_to_trigLPGBTs;
    tId_to_tId       motherboard_to_nTrigLPGBT;

  };


  void OpenGeometry( const std::string& aFilename, const std::string& basePath, Maps& aMaps );

  void OpenRegions( const std::string& aFilename, const std::string& basePath, Maps& aMaps );

  void OpenS1( const std::string& aFilename, const std::string& basePath, Maps& aMaps );

  void OpenChannelAllocation( const std::string& aFilename, const std::string& basePath, Maps& aMaps );

  void OpenBackendMapping( const std::string& aFilename, const std::string& basePath, Maps& aMaps );
}

