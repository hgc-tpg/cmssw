namespace l1thgcmapping {

  inline uint32_t Header(const uint32_t& endcap,
                         const uint32_t& sector,
                         const uint32_t& subsector,
                         const uint32_t& subsystem) {
    if (endcap & ~0x1)
      throw std::runtime_error("Invalid endcap");
    if (sector & ~0x3)
      throw std::runtime_error("Invalid sector");
    if (subsector & ~0x1)
      throw std::runtime_error("Invalid subsector");
    if (subsystem & ~0x3)
      throw std::runtime_error("Invalid subsystem");

    return ((endcap & 0x1) << 31) | ((sector & 0x3) << 29) | ((subsector & 0x1) << 28) | ((subsystem & 0x3) << 26);
  }

  inline uint32_t ObjectType(const uint32_t& object_type) { return ((object_type & 0xF) << 22); }

  inline uint32_t ResetObjectType(const uint32_t& word) { return word & ~(0xF << 22); }

  // ---------------------------------------------------------------------------------------------------------------------------------------------------------

  inline uint32_t get_endcap(const uint32_t& id) { return (id >> 31) & 0x1; }

  inline uint32_t get_sector(const uint32_t& id) { return (id >> 29) & 0x3; }

  inline uint32_t get_subsector(const uint32_t& id) { return (id >> 28) & 0x1; }

  inline uint32_t get_subsystem(const uint32_t& id) { return (id >> 26) & 0x3; }

  inline uint32_t get_object_type(const uint32_t& id) { return (id >> 22) & 0xF; }

  inline std::string get_type_str(const uint32_t& id) {
    uint32_t subsystem(get_subsystem(id)), type(get_object_type(id));
    if (subsystem == 0) {
      if (type == 6)
        return "Scintillator Module";
      if (type == 0)
        return "Silicon Module";
      if (type == 1)
        return "Trigger Cell";
      if (type == 2)
        return "Motherboard";
      if (type == 3)
        return "Trigger lpGBT";
      if (type == 4)
        return "DAQ lpGBT";
      if (type == 5)
        return "Region";
    } else if (subsystem == 1) {
      if (type == 0)
        return "S1";
      if (type == 1)
        return "S1 input channel";
      if (type == 2)
        return "S1 output channel";
      if (type == 3)
        return "Backend fibre";
      if (type == 4)
        return "Backend cluster channel";
      if (type == 5)
        return "Backend tower channel";
    } else if (subsystem == 2) {
      if (type == 0)
        return "S2";
      if (type == 1)
        return "S2 input channel";
      if (type == 2)
        return "S2 output channel";
    }
    throw std::runtime_error("Bad Type");
    return nullptr;
  }

  inline uint32_t add_sector_to_id(const uint32_t& id, const uint32_t& sector) {
    return (id & ~(0x3 << 29)) | ((sector & 0x3) << 29);
  }

  // ---------------------------------------------------------------------------------------------------------------------------------------------------------

  inline uint32_t silicon_module_id(const uint32_t& plane,
                                    const uint32_t& u,
                                    const uint32_t& v,
                                    const uint32_t& endcap = 0,
                                    const uint32_t& sector = 3,
                                    const uint32_t& subsector = 0) {
    if (plane & ~0x3F)
      throw std::runtime_error("Invalid plane");
    if (u & ~0xF)
      throw std::runtime_error("Invalid u");
    if (v & ~0xF)
      throw std::runtime_error("Invalid v");
    return Header(endcap, sector, subsector, 0) | ObjectType(0) | ((plane & 0x3F) << 16) | ((u & 0xF) << 12) |
           ((v & 0xF) << 8);
  }

  inline uint32_t scintillator_module_id(const uint32_t& plane,
                                         const uint32_t& u,
                                         const uint32_t& v,
                                         const uint32_t& endcap = 0,
                                         const uint32_t& sector = 3,
                                         const uint32_t& subsector = 0) {
    if (plane & ~0x3F)
      throw std::runtime_error("Invalid plane");
    if (u & ~0xF)
      throw std::runtime_error("Invalid u");
    if (v & ~0xF)
      throw std::runtime_error("Invalid v");
    return Header(endcap, sector, subsector, 0) | ObjectType(6) | ((plane & 0x3F) << 16) | ((u & 0xF) << 12) |
           ((v & 0xF) << 8);
  }

  inline uint32_t trigger_cell_id(const uint32_t& plane,
                                  const uint32_t& u,
                                  const uint32_t& v,
                                  const uint32_t& cell_u,
                                  const uint32_t& cell_v,
                                  const uint32_t& endcap = 0,
                                  const uint32_t& sector = 3,
                                  const uint32_t& subsector = 0) {
    if (cell_u & ~0xF)
      throw std::runtime_error("Invalid cell_u");
    if (cell_v & ~0xF)
      throw std::runtime_error("Invalid cell_v");
    return silicon_module_id(plane, u, v, endcap, sector, subsector) | ObjectType(1) | ((cell_u & 0xF) << 4) |
           ((cell_v & 0xF) << 0);
  }  // No need to reset object type as Module_id object_type is 0

  inline uint32_t motherboard_id(const uint32_t& plane,
                                 const uint32_t& mbid,
                                 const uint32_t& endcap = 0,
                                 const uint32_t& sector = 3,
                                 const uint32_t& subsector = 0) {
    if (plane & ~0x3F)
      throw std::runtime_error("Invalid plane");
    if (mbid & ~0x1FFF)
      throw std::runtime_error("Invalid mbid");
    return Header(endcap, sector, subsector, 0) | ObjectType(2) | ((plane & 0x3F) << 16) | ((mbid & 0x1FFF) << 3);
  }

  inline uint32_t trigger_lpGBT_id(const uint32_t& plane,
                                   const uint32_t& mbid,
                                   const uint32_t& index,
                                   const uint32_t& endcap = 0,
                                   const uint32_t& sector = 3,
                                   const uint32_t& subsector = 0) {
    if (index & ~0x7)
      throw std::runtime_error("Invalid index");
    return ResetObjectType(motherboard_id(plane, mbid, endcap, sector, subsector)) | ObjectType(3) |
           ((index & 0x7) << 0);
  }

  inline uint32_t daq_lpGBT_id(const uint32_t& plane,
                               const uint32_t& mbid,
                               const uint32_t& index,
                               const uint32_t& endcap = 0,
                               const uint32_t& sector = 3,
                               const uint32_t& subsector = 0) {
    if (index & ~0x7)
      throw std::runtime_error("Invalid index");
    return ResetObjectType(motherboard_id(plane, mbid, endcap, sector, subsector)) | ObjectType(4) |
           ((index & 0x7) << 0);
  }

  inline uint32_t region_id(const uint32_t& plane,
                            const uint32_t& type,
                            const uint32_t& index,
                            const uint32_t& endcap = 0,
                            const uint32_t& sector = 3,
                            const uint32_t& subsector = 0) {
    if (plane & ~0x3F)
      throw std::runtime_error("Invalid plane");
    if (type & ~0x3)
      throw std::runtime_error("Invalid type");
    if (index & ~0x1)
      throw std::runtime_error("Invalid index");
    return Header(endcap, sector, subsector, 0) | ObjectType(5) | ((plane & 0x3F) << 16) | ((type & 0x3) << 1) |
           ((index & 0x1) << 0);
  }

  inline uint32_t patch_id(const uint32_t& is_trigger,
                           const uint32_t& layer,
                           const uint32_t& ID0,
                           const uint32_t& ID1,
                           const uint32_t& ID2,
                           const uint32_t& endcap = 0,
                           const uint32_t& sector = 3,
                           const uint32_t& subsector = 0) {
    if (is_trigger & ~0x1)
      throw std::runtime_error("Invalid trigger/daq flag");
    if (layer & ~0x1)
      throw std::runtime_error("Invalid layer");
    if (ID0 & ~0x7F)
      throw std::runtime_error("Invalid ID0");
    if (ID1 & ~0x7F)
      throw std::runtime_error("Invalid ID1");
    if (ID2 & ~0x3)
      throw std::runtime_error("Invalid ID2");
    return Header(endcap, sector, subsector, 0) | ObjectType(6) | (0x3F << 16) | ((is_trigger & 0x1) << 15) |
           ((layer & 0x1) << 14) | ((ID0 & 0x7F) << 8) | ((ID1 & 0x7F) << 2) | ((ID2 & 0x3) << 0);
  }

  inline uint32_t get_plane(const uint32_t& id) {
    if (get_subsystem(id) != 0)
      throw std::runtime_error("get_planes can only be used in subsystem 0");
    return (id >> 16) & 0x3F;
  }

  inline std::pair<uint32_t, uint32_t> get_module_uv(const uint32_t& id) {
    if (get_subsystem(id) != 0)
      throw std::runtime_error("get_module_uv can only be used in subsystem 0");
    if (!(get_object_type(id) == 0 || get_object_type(id) == 6))
      throw std::runtime_error("get_module_uv can only be used on module-IDs");
    return std::make_pair((id >> 12) & 0xF, (id >> 8) & 0xF);
  }

  inline uint32_t get_motherboard_id(const uint32_t& id) {
    if (get_subsystem(id) != 0)
      throw std::runtime_error("get_motherboard_id can only be used in subsystem 0");
    if (get_object_type(id) < 2 or get_object_type(id) > 3)
      throw std::runtime_error("get_motherboard_id can only be used on motherboard-IDs");
    return (id >> 3) & 0x1FFF;
  }

  inline uint32_t get_region_type(const uint32_t& id) {
    if (get_subsystem(id) != 0)
      throw std::runtime_error("get_region_type can only be used in subsystem 0");
    if (get_object_type(id) < 5)
      throw std::runtime_error("get_region_type can only be used on region-IDs");
    return (id >> 1) & 0x3;
  }

  inline uint32_t get_region_id(const uint32_t& id) {
    if (get_subsystem(id) != 0)
      throw std::runtime_error("get_region_id can only be used in subsystem 0");
    if (get_object_type(id) < 5)
      throw std::runtime_error("get_region_id can only be used on region-IDs");
    return (id >> 0) & 0x1;
  }

  inline uint32_t get_patch_is_trigger(const uint32_t& id) {
    if (get_subsystem(id) != 0)
      throw std::runtime_error("get_patch_is_trigger can only be used in subsystem 0");
    if (get_object_type(id) != 6)
      throw std::runtime_error("get_patch_is_trigger can only be used on patch-IDs");
    return (id >> 15) & 0x1;
  }

  // ---------------------------------------------------------------------------------------------------------------------------------------------------------

  inline uint32_t S1_id(const uint32_t& s1index,
                        const uint32_t& endcap = 0,
                        const uint32_t& sector = 3,
                        const uint32_t& subsector = 0) {
    if (s1index & ~0x3f)
      throw std::runtime_error("Invalid index");
    return Header(endcap, sector, subsector, 1) | ObjectType(0) | ((s1index & 0x3F) << 16);
  }

  inline uint32_t S1_input_channel_id(const uint32_t& s1index,
                                      const uint32_t& index,
                                      const uint32_t& endcap = 0,
                                      const uint32_t& sector = 3,
                                      const uint32_t& subsector = 0) {
    if (index & ~0x7F)
      throw std::runtime_error("Invalid index");
    return S1_id(s1index, endcap, sector, subsector) | ObjectType(1) | ((index & 0x7F) << 0);
  }  // No need to reset object type as S1_id object_type is 0

  inline uint32_t Add_input_channel_to_S1id(const uint32_t& s1id, const uint32_t& index) {
    if (index & ~0x7F)
      throw std::runtime_error("Invalid index");
    return s1id | ObjectType(1) | ((index & 0x7F) << 0);
  }  // No need to reset object type as S1_id object_type is 0

  inline uint32_t S1_output_channel_id(const uint32_t& s1index,
                                       const uint32_t& index,
                                       const uint32_t& endcap = 0,
                                       const uint32_t& sector = 3,
                                       const uint32_t& subsector = 0) {
    if (index & ~0x7F)
      throw std::runtime_error("Invalid index");
    return S1_id(s1index, endcap, sector, subsector) | ObjectType(2) | ((index & 0x7F) << 0);
  }  // No need to reset object type as S1_id object_type is 0

  inline uint32_t Add_output_channel_to_S1id(const uint32_t& s1id, const uint32_t& index) {
    if (index & ~0x7F)
      throw std::runtime_error("Invalid index");
    return s1id | ObjectType(2) | ((index & 0x7F) << 0);
  }  // No need to reset object type as S1_id object_type is 0

  inline uint32_t S1_backend_fibre_id(const uint32_t& s1index,
                                      const uint32_t& tmindex,
                                      const uint32_t& fibre,
                                      const uint32_t& endcap = 0,
                                      const uint32_t& sector = 3,
                                      const uint32_t& subsector = 0) {
    if (tmindex & ~0x1F)
      throw std::runtime_error("Invalid tmindex");
    if (fibre & ~0x7)
      throw std::runtime_error("Invalid fibre");
    return S1_id(s1index, endcap, sector, subsector) | ObjectType(3) | ((tmindex & 0x1F) << 11) | ((fibre & 0x7) << 8);
  }  // No need to reset object type as S1_id object_type is 0

  inline uint32_t Add_backend_fibre_to_S1id(const uint32_t& s1id, const uint32_t& tmindex, const uint32_t& fibre) {
    if (tmindex & ~0x1F)
      throw std::runtime_error("Invalid tmindex");
    if (fibre & ~0x7)
      throw std::runtime_error("Invalid fibre");
    return s1id | ObjectType(3) | ((tmindex & 0x1F) << 11) | ((fibre & 0x7) << 8);
  }  // No need to reset object type as S1_id object_type is 0

  inline uint32_t S1_backend_cluster_channel_id(const uint32_t& s1index,
                                                const uint32_t& tmindex,
                                                const uint32_t& fibre,
                                                const uint32_t& channel,
                                                const uint32_t& endcap = 0,
                                                const uint32_t& sector = 3,
                                                const uint32_t& subsector = 0) {
    if (channel & ~0x3)
      throw std::runtime_error("Invalid channel");
    return ResetObjectType(S1_backend_fibre_id(s1index, tmindex, fibre, endcap, sector, subsector)) | ObjectType(4) |
           ((channel & 0x3) << 6);
  }

  inline uint32_t Add_backend_cluster_channel_to_S1id(const uint32_t& s1id,
                                                      const uint32_t& tmindex,
                                                      const uint32_t& fibre,
                                                      const uint32_t& channel) {
    if (fibre & ~0x7)
      throw std::runtime_error("Invalid fibre");
    if (channel & ~0x3)
      throw std::runtime_error("Invalid channel");
    return s1id | ObjectType(4) | ((tmindex & 0x1F) << 11) | ((fibre & 0x7) << 8) | ((channel & 0x3) << 6);
  }  // No need to reset object type as S1_id object_type is 0

  // uint32_t S1_backend_cluster_frame_id( s1index , tmindex , fibre , channel , frame , endcap = 0 , sector = 3 , subsector = 0 ){
  // return ResetObjectType( S1_backend_cluster_channel_id( s1index , tmindex , fibre , channel , endcap , sector , subsector ) ) | ObjectType( 5 ) | ( ( frame & 0x7F ) << 0 ); }

  // uint32_t Add_backend_cluster_frame_to_S1id( s1id , tmindex , fibre , channel , frame ){
  // return s1id | ObjectType( 5 ) | ( ( tmindex & 0x1F ) << 11 ) | ( ( fibre & 0x3 ) << 9 ) | ( ( channel & 0x3 ) << 7 ) | ( ( frame & 0x7F ) << 0 ); } // No need to reset object type as S1_id object_type is 0

  // uint32_t S1_backend_tower_channel_id( s1index , tmindex , fibre , channel , endcap = 0 , sector = 3 , subsector = 0 ){
  // return ResetObjectType( S1_backend_fibre_id( s1index , tmindex , fibre , endcap , sector , subsector ) ) | ObjectType( 5 ) | ( ( channel & 0x3 ) << 6 ); }

  // uint32_t S1_backend_tower_frame_id( s1index , tmindex , fibre , channel , frame , endcap = 0 , sector = 3 , subsector = 0 ){
  // return ResetObjectType( S1_backend_tower_channel_id( s1index , tmindex , fibre , channel , endcap , sector , subsector ) ) | ObjectType( 7 ) | ( ( frame & 0x7F ) << 0 ); }

  inline uint32_t get_s1index(const uint32_t& id) {
    if (get_subsystem(id) != 1)
      throw std::runtime_error("get_s1index can only be used in subsystem 1");
    return (id >> 16) & 0x3F;
  }

  inline uint32_t get_backend_fibre(const uint32_t& id) {
    if (get_subsystem(id) != 1)
      throw std::runtime_error("get_backend_fibre can only be used in subsystem 1");
    if (get_object_type(id) < 3)
      throw std::runtime_error("get_backend_fibre can only be used on region-IDs");
    return (id >> 8) & 0x7;
  }

  inline uint32_t get_backend_channel(const uint32_t& id) {
    if (get_subsystem(id) != 1)
      throw std::runtime_error("get_backend_channel can only be used in subsystem 1");
    if (get_object_type(id) < 4)
      throw std::runtime_error("get_backend_channel can only be used on region-IDs");
    return (id >> 6) & 0x3;
  }

  // ---------------------------------------------------------------------------------------------------------------------------------------------------------

  inline uint32_t S2_id(const uint32_t& tmindex, const uint32_t& endcap = 0, const uint32_t& sector = 3) {
    if (tmindex & ~0x1F)
      throw std::runtime_error("Invalid tmindex");
    return Header(endcap, sector, 0, 2) | ObjectType(0) | ((tmindex & 0x1F) << 11);
  }

  inline uint32_t S2_input_channel_id(const uint32_t& tmindex,
                                      const uint32_t& index,
                                      const uint32_t& endcap = 0,
                                      const uint32_t& sector = 3) {
    if (index & ~0x7F)
      throw std::runtime_error("Invalid index");
    return S2_id(tmindex, endcap, sector) | ObjectType(1) | ((index & 0x7F) << 0);
  }  // No need to reset object type as S2_id object_type is 0

  inline uint32_t Add_input_channel_to_S2id(const uint32_t& s2id, const uint32_t& index) {
    if (index & ~0x7F)
      throw std::runtime_error("Invalid index");
    return s2id | ObjectType(1) | ((index & 0x7F) << 0);
  }  // No need to reset object type as S2_id object_type is 0

  inline uint32_t S2_output_channel_id(const uint32_t& tmindex,
                                       const uint32_t& index,
                                       const uint32_t& endcap = 0,
                                       const uint32_t& sector = 3) {
    if (index & ~0x7F)
      throw std::runtime_error("Invalid index");
    return S2_id(tmindex, endcap, sector) | ObjectType(2) | ((index & 0x7F) << 0);
  }  // No need to reset object type as S2_id object_type is 0

  inline uint32_t Add_output_channel_to_S2id(const uint32_t& s2id, const uint32_t& index) {
    if (index & ~0x7F)
      throw std::runtime_error("Invalid index");
    return s2id | ObjectType(2) | ((index & 0x7F) << 0);
  }  // No need to reset object type as S2_id object_type is 0

  inline uint32_t get_s2tm(const uint32_t& id) {
    if (get_subsystem(id) != 2)
      throw std::runtime_error("get_s2tm can only be used in subsystem 2");
    return (id >> 11) & 0x1F;
  }

  inline uint32_t get_io_index(const uint32_t& id) {
    if (get_subsystem(id) == 0)
      throw std::runtime_error("get_io_index can only be used on S1 and S2 subsystems");
    if (get_object_type(id) != 1 and get_object_type(id) != 2)
      throw std::runtime_error("get_io_index can only be used on S1 or S2 I/O IDs");
    return (id >> 0) & 0x7F;
  }
  // ---------------------------------------------------------------------------------------------------------------------------------------------------------
}  // namespace l1thgcmapping
