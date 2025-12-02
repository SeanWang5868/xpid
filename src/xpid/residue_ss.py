"""
ss.py
Secondary structure detection and indexing.
Encapsulates logic to map residues to secondary structure types and unique region IDs.
"""
import gemmi
from typing import Dict, List, Tuple

def build_index(structure: gemmi.Structure) -> Dict[str, List[Tuple[int, int, str, int]]]:
    """
    Builds an index of secondary structures with unique region IDs.
    
    Logic:
    - Helices: Classified into G (3_10), I (Pi), H (Alpha/Others).
    - Sheets: Classified as E (Extended).
    - UIDs: A unique integer ID is assigned to each distinct SS element.
    
    Returns: 
        { 'ChainName': [ (StartSeqNum, EndSeqNum, 'TypeChar', RegionUID), ... ] }
    """
    ss_index = {}
    region_uid_counter = 1

    def add_range(chain_name, start_num, end_num, ss_type, uid):
        if chain_name not in ss_index:
            ss_index[chain_name] = []
        ss_index[chain_name].append((start_num, end_num, ss_type, uid))

    # 1. Process Helices
    if hasattr(structure, 'helices'):
        for helix in structure.helices:
            try:
                # Use seqid.num (int) for robust comparison
                start_num = helix.start.res_id.seqid.num
                end_num = helix.end.res_id.seqid.num
                chain = helix.start.chain_name
                
                # Determine Helix Type
                h_class = helix.pdb_helix_class
                ss_code = 'H' # Default Alpha
                
                if h_class == 5:
                    ss_code = 'G' # 3_10 Helix
                elif h_class == 3:
                    ss_code = 'I' # Pi Helix
                
                add_range(chain, start_num, end_num, ss_code, region_uid_counter)
                region_uid_counter += 1
            except Exception:
                continue

    # 2. Process Sheets (Strands)
    if hasattr(structure, 'sheets'):
        for sheet in structure.sheets:
            for strand in sheet.strands:
                try:
                    start_num = strand.start.res_id.seqid.num
                    end_num = strand.end.res_id.seqid.num
                    chain = strand.start.chain_name
                    
                    add_range(chain, start_num, end_num, 'E', region_uid_counter)
                    region_uid_counter += 1
                except Exception:
                    continue

    return ss_index

def get_info(chain_name: str, res_seq_num: int, ss_index: Dict) -> Tuple[str, int]:
    """
    Queries secondary structure info for a specific residue.
    
    Returns:
        (TypeChar, RegionUID). 
        If no structure found (Coil), returns ('C', -1).
    """
    if chain_name not in ss_index:
        return ('C', -1)
    
    for start_num, end_num, ss_type, uid in ss_index[chain_name]:
        if start_num <= res_seq_num <= end_num:
            return (ss_type, uid)
            
    return ('C', -1)