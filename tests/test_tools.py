import pytest
from ic4popsyn import tools, constants

pytestmark = pytest.mark.skipif(constants.G != 3.920659e8, reason='Different G constanst is set up')
def test_a2p():
    #solar system: given the separation get the period
    years, days = tools.a2p(214.95, 1., 3.32946e-6)
    assert abs(round(years,1)-1) < 1e-10 and abs(round(days,2)-365.25) < 1e-10

def test_p2a():
    #solar system: given the period get the separation
    Rsun = tools.p2a(365.25, 1., 3.33e-6)
    assert abs(round(Rsun,2)-214.95) < 1e-10
