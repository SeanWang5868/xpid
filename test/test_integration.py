import gemmi
from xpid import core, config

def test_core_detection_empty():
    st = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    st.add_model(model)
    model.add_chain(chain)
    
    # Empty structure should return empty list
    res = core.detect_interactions_in_structure(st, "test", {}, model_mode=0)
    assert len(res) == 0

def test_model_selection():
    st = gemmi.Structure()
    m1 = gemmi.Model("1")
    m2 = gemmi.Model("2")
    st.add_model(m1)
    st.add_model(m2)
    
    # Dummy atoms to prevent immediate return
    c1 = gemmi.Chain("A")
    m1.add_chain(c1)
    
    # Test valid index
    core.detect_interactions_in_structure(st, "test", {}, model_mode=0)
    # Test 'all'
    core.detect_interactions_in_structure(st, "test", {}, model_mode='all')