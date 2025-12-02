import gemmi
import os

def load_and_prep_structure(filepath):
    """读取结构并进行预处理（如果需要加氢）"""
    try:
        structure = gemmi.read_structure(str(filepath))
        
        # 检查是否有氢原子
        has_hydrogen = False
        for model in structure:
            for chain in model:
                for res in chain:
                    for atom in res:
                        if atom.element.name == 'H':
                            has_hydrogen = True
                            break
                    if has_hydrogen: break
                if has_hydrogen: break
        
        # 如果没有氢原子，尝试使用 gemmi 添加 (简单的几何加氢)
        if not has_hydrogen:
            # print(f"Info: No hydrogens found in {filepath}, attempting to add...")
            structure.setup_entities()
            structure.assign_label_seq_id()
            # 注意：Gemmi 的加氢相对基础，如果需要精确加氢建议外部使用 reduce
            structure.add_hydrogens()
            
        return structure
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None