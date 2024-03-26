# tests.py performs basic unit testing for the optimisation module
# command is pytest tests.py

from optim_module import OptimModule

optim_module = OptimModule("test_optimisation.json")


def cavity_optimisation():

    thicknesses = optim_module.perform_optimisation("minimize")
    return thicknesses


def test():
    assert optim_module.check_targets(cavity_optimisation())
