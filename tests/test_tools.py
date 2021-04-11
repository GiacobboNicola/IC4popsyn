import numpy as np
from ic4popsyn import tools

def test_a2p():
    #solar system: given the separation get the period
    years, days = tools.a2p(214.95, 1, 3e-6)
    assert abs(round(years,1)-1) < 1e-10 and abs(round(days,2)-365.24) < 1e-10

def test_p2a():
    #solar system: given the period get the separation
    Rsun = tools.p2a(365.24, 1, 3e-6)
    assert abs(round(Rsun,2)-214.95) < 1e-10
